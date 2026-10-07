#!/usr/bin/env python3
"""Cross-check the §4 evaluator against CAFA-evaluator (``cafaeval``) — TODO P0.

Our protein-centric metrics (``temporal_benchmark`` + ``resampling``) are home
grown. This rebuilds the committed ``ipr_manual`` ablation cell with
``ablation.py``'s own functions — same cohort, same truth, same p-score
transfer, same IC floor — exports each rung to CAFA format and scores it with
the reference implementation of Piovesan et al. 2024 (Bioinformatics Advances,
doi:10.1093/bioadv/vbae043), then compares, per rung x aspect x IC floor.

What is handed to cafaeval, and what it recomputes itself
-----------------------------------------------------------
* **Truth.** The *unpropagated* t1 experimental terms of each cohort protein
  in the aspect. cafaeval propagates them over its own parse of the same
  ``go-basic.obo`` (``is_a`` + ``part_of``, obsolete terms dropped), so our
  truth propagation is checked rather than assumed. The cohort itself — CAFA
  no-knowledge targets with a domain — is ours: cafaeval has no notion of it.
* **Predictions.** Our transferred per-protein scores, already propagated
  (the p-score sums over ancestors). cafaeval max-propagates them again; for a
  consistent prediction that is a no-op, which the term-level comparison
  below verifies.
* **IC weights.** Our marginal t0 IC, passed as cafaeval's information
  accretion file, so its IA-weighted ``ru``/``mi`` are our S ingredients.
* **IC floor.** cafaeval has no floor. Its ``-no_orphans`` option works by
  restricting ``Graph.toi`` (the evaluated term indices); the floor is applied
  through that same attribute — ``toi`` is set to the terms ``filter_by_ic``
  keeps, aspect roots excluded — and the ground-truth file holds only the
  floor's cohort, because cafaeval counts every ground-truth protein. That is
  exactly ``ablation.floor_cell``: both sides restricted, emptied proteins
  leave the cohort.

Two comparisons
---------------
1. **Same thresholds** (the implementation check). cafaeval's
   ``evaluate_prediction`` takes any threshold array, so it is run at the
   cell's own sweep (``_candidate_thresholds``) and at the 0.001 grid, and the
   precision / recall / coverage / S curves are compared point by point with
   ours at those thresholds (``curve_max_abs_diff``). One point is left out:
   ``tau = 0``. In CAFA format a zero score *is* "not predicted" and
   ``pred >= 0`` would predict every term of the ontology, so our
   predict-everything-scored point cannot be expressed.
2. **As reported.** Our committed numbers against cafaeval run the way
   ``cafa_eval`` runs: thresholds ``arange(th_step, 1, th_step)`` with
   ``th_step = 0.001`` (ten times its default, so the grid is not the
   bottleneck), rows with zero coverage dropped, F_max = first max of ``f``.
   S_min is the minimum of its IA-weighted ``s_w`` column: cafaeval 1.3.0's own
   ``evaluation_best_s`` minimises the *unweighted* (term-count) S even when an
   IA file is given. cafaeval reports no area, so AUPRC is our trapezoid
   (``temporal_benchmark.auprc``) over its curve.

Regeneration::

    uv run python validation/check_evaluator_cafaeval.py \\
        --t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz \\
        --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \\
        --run-dir results/ablation-replacedby/ipr_manual \\
        --work-dir results/ablation-cafaeval/ipr_manual
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent


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
ASPECT_NAMESPACE = {
    aspect: namespace for namespace, aspect in tb.NAMESPACE_TO_ASPECT.items()
}


# --------------------------------------------------------------------------- #
# CAFA-format export                                                           #
# --------------------------------------------------------------------------- #
def write_ground_truth(path: Path, truth: Mapping[str, Iterable[str]]) -> int:
    """Write ``protein<TAB>term`` lines (sorted); returns the number written."""
    n = 0
    with path.open("w") as handle:
        for protein in sorted(truth):
            for term in sorted(truth[protein]):
                handle.write(f"{protein}\t{term}\n")
                n += 1
    return n


def write_predictions(path: Path, preds: Mapping[str, Mapping[str, float]]) -> int:
    """Write ``protein<TAB>term<TAB>score`` lines (sorted); returns the number written.

    CAFA scores live in ``(0, 1]``: cafaeval stores predictions in a zero-filled
    matrix and its thresholds start above zero, so a score of 0 is "not
    predicted" and is not written. Scores are written with ``repr`` precision so
    that every ``score >= tau`` comparison comes out as in our evaluator.
    """
    n = 0
    with path.open("w") as handle:
        for protein in sorted(preds):
            for term, score in sorted(preds[protein].items()):
                if not 0.0 <= score <= 1.0:
                    raise ValueError(
                        f"CAFA scores must lie in [0, 1]: {protein} {term} {score!r}"
                    )
                if score > 0.0:
                    handle.write(f"{protein}\t{term}\t{float(score)!r}\n")
                    n += 1
    return n


def write_information_accretion(path: Path, ic: Mapping[str, float]) -> int:
    """Write ``term<TAB>weight`` lines — our IC as cafaeval's IA weights."""
    with path.open("w") as handle:
        for term in sorted(ic):
            handle.write(f"{term}\t{float(ic[term])!r}\n")
    return len(ic)


def raw_truth(
    t1_exp_map: Mapping[str, set[str]],
    cohort: Iterable[str],
    aspect: str,
    term_aspect: Mapping[str, str],
) -> dict[str, set[str]]:
    """The cohort's unpropagated t1 terms in one aspect, for cafaeval to propagate."""
    out: dict[str, set[str]] = {}
    for protein in cohort:
        terms = {t for t in t1_exp_map.get(protein, ()) if term_aspect.get(t) == aspect}
        if terms:
            out[protein] = terms
    return out


def floor_terms(
    terms: Iterable[str], ic: Mapping[str, float], min_ic: float
) -> set[str]:
    """The terms an IC floor keeps — ``filter_by_ic``'s own predicate — minus roots."""
    kept = tb.filter_by_ic({"": set(terms)}, ic, min_ic).get("", set())
    return kept - tb.ASPECT_ROOTS


# --------------------------------------------------------------------------- #
# Metrics from a curve                                                         #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Curve:
    """A protein-centric curve: one value per threshold, thresholds ascending."""

    tau: np.ndarray
    precision: np.ndarray
    recall: np.ndarray
    coverage: np.ndarray
    s: np.ndarray


def curve_summary(curve: Curve) -> dict[str, float]:
    """F_max (first max over ascending tau), its tau and coverage, S_min, AUPRC.

    The same conventions ``resampling.panel_metrics`` uses, so that applied to
    our own curve it returns our numbers and applied to cafaeval's it returns
    the numbers cafaeval's curve implies.
    """
    if len(curve.tau) == 0:
        return {
            "f_max": 0.0,
            "f_max_tau": 0.0,
            "coverage_at_fmax": 0.0,
            "s_min": math.nan,
            "s_min_tau": math.nan,
            "auprc": 0.0,
        }
    p, r = curve.precision, curve.recall
    denom = p + r
    f = np.where(denom > 0, 2 * p * r / np.maximum(denom, 1e-300), 0.0)
    best = int(np.argmax(f))
    s_best = int(np.argmin(curve.s))
    return {
        "f_max": float(f[best]),
        "f_max_tau": float(curve.tau[best]) if f[best] > 0 else 0.0,
        "coverage_at_fmax": float(curve.coverage[best]),
        "s_min": float(curve.s[s_best]),
        "s_min_tau": float(curve.tau[s_best]),
        "auprc": tb.auprc(list(zip(curve.tau, p, r))),
    }


def max_abs_diff(a: Curve, b: Curve) -> float:
    """Largest pointwise disagreement between two curves on the same thresholds."""
    if not np.array_equal(a.tau, b.tau):
        raise ValueError("curves are on different thresholds")
    if len(a.tau) == 0:
        return 0.0
    return float(
        max(
            np.max(np.abs(a.precision - b.precision)),
            np.max(np.abs(a.recall - b.recall)),
            np.max(np.abs(a.coverage - b.coverage)),
            np.max(np.abs(a.s - b.s)),
        )
    )


def our_curve(pred_f, truth_f, ic, taus) -> Curve:
    """Our evaluator's curve (``resampling.panel_curve``) at the given thresholds."""
    panel = rs.build_panel(pred_f, truth_f, ic, taus)
    c = rs.panel_curve(panel)
    return Curve(panel.thresholds, c.precision, c.recall, c.coverage, c.s)


def cafaeval_curve(df) -> Curve:
    """cafaeval's ``evaluate_prediction`` frame (one namespace) as a :class:`Curve`."""
    df = df.sort_values("tau")
    return Curve(
        df["tau"].to_numpy(dtype=float),
        df["pr"].to_numpy(dtype=float),
        df["rc"].to_numpy(dtype=float),
        df["cov"].to_numpy(dtype=float),
        df["s_w"].to_numpy(dtype=float),
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def manifest_hash_mismatches(
    manifest_dir: Path, expected: Mapping[str, Path]
) -> list[str]:
    """Inputs whose sha256 differs from what the training runs recorded.

    ``expected`` maps a manifest input ``role`` (``gaf``, ``go_obo``,
    ``domain_annotations``) to the file this check reads for it. A role no
    manifest records is a problem too — a wrong or empty ``manifest_dir``
    would otherwise pass without comparing anything.
    """
    hashes = {role: sha256(path) for role, path in expected.items()}
    problems = []
    seen: set[str] = set()
    for manifest in sorted(manifest_dir.glob("*.json")):
        for entry in json.loads(manifest.read_text()).get("inputs", []):
            want = hashes.get(entry.get("role"))
            if want is None:
                continue
            seen.add(entry["role"])
            if entry.get("sha256") != want:
                problems.append(
                    f"{manifest.name}: {entry['role']} {entry['sha256'][:12]} "
                    f"!= {expected[entry['role']]} {want[:12]}"
                )
    for role in sorted(set(expected) - seen):
        problems.append(f"no manifest under {manifest_dir} records role {role}")
    return problems


# --------------------------------------------------------------------------- #
# The cafaeval side (one worker per method x floor)                            #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Job:
    """Shared state the forked workers inherit.

    ``ontologies`` and ``gts`` (one parse per IC floor) are parsed once, before
    the fork; a worker only overwrites ``Graph.toi`` / ``toi_ia``, and does so
    for every aspect before reading them.
    """

    ontologies: dict
    gts: dict[float, dict]
    pred_files: dict[str, Path]
    cells: dict[tuple[str, str, float], tuple[dict, dict, list[float]]]
    ic: dict[str, float]
    grid: np.ndarray
    n_threads: int


_JOB: Job | None = None


def parse_cafa(
    obo: Path, ia_file: Path, gt_files: Mapping[float, Path]
) -> tuple[dict, dict[float, dict]]:
    """cafaeval's own parse of the ontology (with IA) and of each ground truth.

    ``is_a`` + ``part_of``, orphans (roots) kept — the evaluated terms are set
    per aspect and floor in :func:`_score_cell`. Truth is propagated here.
    """
    from cafaeval.parser import gt_parser, obo_parser

    ontologies = obo_parser(str(obo), ("is_a", "part_of"), str(ia_file), True)
    gts = {k: gt_parser(str(path), ontologies) for k, path in gt_files.items()}
    return ontologies, gts


def _score_cell(task: tuple[str, float]) -> list[dict]:
    """Score one method at one IC floor with cafaeval, every aspect."""
    from cafaeval.evaluation import evaluate_prediction
    from cafaeval.parser import pred_parser
    from loguru import logger

    name, min_ic = task
    job = _JOB
    assert job is not None
    ontologies, gts = job.ontologies, job.gts[min_ic]
    preds = pred_parser(str(job.pred_files[name]), ontologies, gts, "max")
    rows = []
    for aspect in ASPECTS:
        cell = job.cells.get((name, aspect, min_ic))
        if cell is None:
            continue
        truth_f, pred_f, taus = cell
        ont = ontologies[ASPECT_NAMESPACE[aspect]]
        gt = gts[ASPECT_NAMESPACE[aspect]]
        index = {t: v["index"] for t, v in ont.terms_dict.items()}
        toi_terms = floor_terms(index, job.ic, min_ic)
        ont.toi = np.array(sorted(index[t] for t in toi_terms), dtype=int)
        ont.toi_ia = ont.toi[ont.ia[ont.toi] > 0]
        terms = [ont.terms_list[i]["id"] for i in ont.toi]

        # Term-level agreement: cafaeval's propagated truth and predictions
        # on the evaluated terms against ours.
        truth_mismatch = 0
        for protein, row in gt.ids.items():
            theirs = {terms[j] for j in np.nonzero(gt.matrix[row, ont.toi])[0]}
            truth_mismatch += len(theirs ^ truth_f.get(protein, set()))
        truth_mismatch += len(set(truth_f) - set(gt.ids))
        pred_mismatch = 0
        prediction = preds.get(ASPECT_NAMESPACE[aspect])
        for protein, row in gt.ids.items():
            ours = {t: s for t, s in pred_f.get(protein, {}).items() if s > 0}
            theirs = {}
            if prediction is not None and protein in prediction.ids:
                values = prediction.matrix[prediction.ids[protein], ont.toi]
                theirs = {terms[j]: float(values[j]) for j in np.nonzero(values)[0]}
            pred_mismatch += sum(
                1 for t in ours.keys() | theirs.keys() if ours.get(t) != theirs.get(t)
            )

        ns = ASPECT_NAMESPACE[aspect]
        same_taus = np.asarray([t for t in taus if t > 0], dtype=float)
        if prediction is None or not len(same_taus):
            raise RuntimeError(f"{name} {aspect}: nothing for cafaeval to score")
        cafa_same = cafaeval_curve(
            evaluate_prediction(
                {ns: prediction}, gts, ontologies, same_taus, "cafa", job.n_threads
            )
        )
        cafa_grid_df = evaluate_prediction(
            {ns: prediction}, gts, ontologies, job.grid, "cafa", job.n_threads
        )
        cafa_grid = cafaeval_curve(cafa_grid_df)
        diff = max(
            max_abs_diff(our_curve(pred_f, truth_f, job.ic, same_taus), cafa_same),
            max_abs_diff(our_curve(pred_f, truth_f, job.ic, job.grid), cafa_grid),
        )
        # As cafa_eval reports it: zero-coverage thresholds dropped.
        reported = cafaeval_curve(cafa_grid_df[cafa_grid_df["cov"] > 0])
        rows.append(
            {
                "aspect": aspect,
                "min_ic": min_ic,
                "method": name,
                "n_eval_cafaeval": len(gt.ids),
                "truth_mismatch": truth_mismatch,
                "prediction_mismatch": pred_mismatch,
                "curve_max_abs_diff": diff,
                **{f"{k}_cafaeval": v for k, v in curve_summary(reported).items()},
            }
        )
        logger.info(
            f"  [{aspect} IC>={min_ic:g}] {name}: curve |Δ|<={diff:.1e}, "
            f"truth/prediction mismatches {truth_mismatch}/{pred_mismatch}"
        )
    return rows


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
        description="Cross-check the §4 ablation evaluator against CAFA-evaluator."
    )
    parser.add_argument("--t0-gaf", type=Path, required=True)
    parser.add_argument("--t1-gaf", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--work-dir",
        type=Path,
        required=True,
        help="Where the CAFA-format files go (large; keep it out of git)",
    )
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
        "--manifest-dir", type=Path, default=_HERE / "ablation_manifests"
    )
    parser.add_argument(
        "--committed-metrics", type=Path, default=_HERE / "ablation_metrics.tsv"
    )
    parser.add_argument(
        "--output", type=Path, default=_HERE / "evaluator_crosscheck.tsv"
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-9,
        help="Sanity gate: max |recomputed - committed| for any metric. Not zero: "
        "IC sums run over sets, so their last bit depends on the hash seed",
    )
    parser.add_argument("--th-step", type=float, default=0.001)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--threads", type=int, default=2, help="cafaeval threads per worker"
    )
    args = parser.parse_args()
    ic_floors = [0.0, 2.0, 4.0]

    # The t1 GAF is not a training input, so no manifest records it; the
    # committed-metrics gate below is what catches a different one.
    problems = manifest_hash_mismatches(
        args.manifest_dir,
        {
            "gaf": args.t0_gaf,
            "go_obo": args.go_ontology,
            "domain_annotations": args.interpro,
        },
    )
    problems += abl.run_setting_mismatches(
        args.run_dir,
        abl.LADDER,
        {"domain_key": args.domain_key, "evidence_filter": args.evidence_filter},
    )
    if problems:
        for problem in problems:
            logger.error(f"Input mismatch — {problem}")
        return 1
    logger.info(f"✓ t0 GAF, GO ontology and domain file match {args.manifest_dir}")

    # --------------------------------------------- our side, as ablation.py --
    inputs = abl.load_benchmark_inputs(
        args.t0_gaf,
        args.t1_gaf,
        args.interpro,
        args.go_ontology,
        args.domain_key,
        args.evidence_filter,
    )
    methods, _, _ = abl.build_methods(args.run_dir, inputs, "pscore", args.domain_key)
    committed = pd.read_csv(args.committed_metrics, sep="\t")
    committed_lookup = {
        (r["aspect"], float(r["min_ic"]), r["method"]): r
        for _, r in committed.iterrows()
    }

    cells: dict[tuple[str, str, float], tuple[dict, dict, list[float]]] = {}
    ours_rows: dict[tuple[str, float, str], dict] = {}
    for name, preds in methods.items():
        for aspect in ASPECTS:
            pred_a = tb.restrict_to_aspect(preds, aspect, inputs.term_aspect)
            for min_ic in ic_floors:
                cell = abl.floor_cell(
                    inputs.benchmark[aspect], pred_a, inputs.ic, min_ic
                )
                truth_f, pred_f, taus = cell
                if not truth_f:
                    continue
                observed = rs.panel_metrics(
                    rs.build_panel(pred_f, truth_f, inputs.ic, taus)
                )
                com = committed_lookup[(aspect, min_ic, name)]
                for key in ("f_max", "f_max_tau", "s_min", "auprc", "coverage_at_fmax"):
                    if abs(observed[key] - com[key]) > args.tolerance:
                        logger.error(
                            f"[{aspect} IC>={min_ic:g}] {name}: recomputed {key} "
                            f"{observed[key]!r} != committed {com[key]!r}"
                        )
                        return 1
                cells[(name, aspect, min_ic)] = cell
                ours_rows[(aspect, min_ic, name)] = observed
    logger.info(f"✓ Recomputed {len(cells)} cells reproduce {args.committed_metrics}")

    # ------------------------------------------------------- CAFA export --
    args.work_dir.mkdir(parents=True, exist_ok=True)
    ia_file = args.work_dir / "information_accretion.tsv"
    write_information_accretion(ia_file, inputs.ic)
    gt_files = {}
    for min_ic in ic_floors:
        truth = {}
        for aspect in ASPECTS:
            cohort = tb.filter_by_ic(inputs.benchmark[aspect], inputs.ic, min_ic)
            for protein, terms in raw_truth(
                inputs.t1_exp_map, cohort, aspect, inputs.term_aspect
            ).items():
                truth.setdefault(protein, set()).update(terms)
        gt_files[min_ic] = args.work_dir / f"ground_truth_ic{min_ic:g}.tsv"
        n = write_ground_truth(gt_files[min_ic], truth)
        logger.info(f"  ground truth IC>={min_ic:g}: {n:,} raw annotations")
    pred_dir = args.work_dir / "predictions"
    pred_dir.mkdir(exist_ok=True)
    pred_files = {}
    for name, preds in methods.items():
        pred_files[name] = pred_dir / f"{name}.tsv"
        n = write_predictions(pred_files[name], preds)
        n_zero = sum(1 for p in preds.values() for s in p.values() if s == 0.0)
        logger.info(
            f"  {name}: {n:,} predictions written, {n_zero:,} zero scores left out"
        )

    # ----------------------------------------------------------- cafaeval --
    logger.info("cafaeval: parsing the ontology and propagating the truth...")
    ontologies, gts = parse_cafa(args.go_ontology, ia_file, gt_files)
    _JOB = Job(
        ontologies=ontologies,
        gts=gts,
        pred_files=pred_files,
        cells=cells,
        ic=inputs.ic,
        grid=np.arange(args.th_step, 1, args.th_step),
        n_threads=args.threads,
    )
    logger.info(f"Scoring {len(methods)} methods with cafaeval...")
    tasks = [(name, min_ic) for name in methods for min_ic in ic_floors]
    with mp.get_context("fork").Pool(args.workers) as pool:
        theirs = {
            (row["aspect"], row["min_ic"], row["method"]): row
            for rows in pool.imap_unordered(_score_cell, tasks)
            for row in rows
        }

    out = []
    for aspect in ASPECTS:
        for min_ic in ic_floors:
            for name in methods:
                key = (aspect, min_ic, name)
                if key not in ours_rows:
                    continue
                ours, row = ours_rows[key], theirs[key]
                merged = {
                    "aspect": aspect,
                    "min_ic": min_ic,
                    "method": name,
                    "n_eval_ours": ours["n_eval_proteins"],
                    "n_eval_cafaeval": row["n_eval_cafaeval"],
                    "truth_mismatch": row["truth_mismatch"],
                    "prediction_mismatch": row["prediction_mismatch"],
                    "curve_max_abs_diff": row["curve_max_abs_diff"],
                }
                for metric in (
                    "f_max",
                    "f_max_tau",
                    "coverage_at_fmax",
                    "s_min",
                    "s_min_tau",
                    "auprc",
                ):
                    merged[f"{metric}_ours"] = ours[metric]
                    merged[f"{metric}_cafaeval"] = row[f"{metric}_cafaeval"]
                    if not metric.endswith("_tau"):
                        merged[f"{metric}_diff"] = (
                            row[f"{metric}_cafaeval"] - ours[metric]
                        )
                out.append(merged)
    df = pd.DataFrame(out)
    df.to_csv(args.output, sep="\t", index=False)
    logger.info(f"✓ Cross-check ({len(df)} cells): {args.output}")
    for col in (
        "truth_mismatch",
        "prediction_mismatch",
        "curve_max_abs_diff",
        "f_max_diff",
        "coverage_at_fmax_diff",
        "s_min_diff",
        "auprc_diff",
    ):
        logger.info(f"  max |{col}| = {df[col].abs().max():.3g}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
