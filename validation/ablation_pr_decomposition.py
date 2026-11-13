#!/usr/bin/env python3
"""Decompose each ablation rung's F_max into precision (PPV) and recall.

For the relative-inference report: the committed ablation table
(``validation/ablation_metrics.tsv``, VALIDATION_PLAN §4) shows relative
inference lowering F_max, but F_max alone cannot say *how* — a precision gain
swamped by a recall loss looks the same as losing both. This script rebuilds
the §4 evaluation exactly (same inputs, same cohort, same panels; no bootstraps,
no permutations) and, for every rung × aspect × IC-floor cell, reports the CAFA
precision and recall vectors' values **at the F_max operating point** — PPV
averaged over only the proteins with ≥1 prediction at that cutoff, sensitivity
averaged over the whole cohort, exactly as ``resampling.panel_metrics``
maximises F.

Two fixed-operating-point columns make rungs comparable off their own F_max:
``ppv_at_ref_sensitivity`` is the rung's PPV, read off its own curve by linear
interpolation, at the sensitivity the reference rung (``--reference-rung``,
default ``supra_input_relative``) reaches at *its* F_max; ``f_at_budget`` is the
rung's F at the number of (protein, term) predictions the reference rung makes
at its F_max. Both are NaN when the rung's curve never reaches that point.

Hard sanity gate: the recomputed ``f_max``/``auprc``/``n_eval_proteins`` must
match the committed ``ablation_metrics.tsv`` to within ``--tolerance`` in every
cell, or nothing is written. ``--min-tau``/``--anchor-origin`` select the
scoring convention; the defaults are ``resampling.panel_metrics``'s CAFA
conventions (F_max at tau >= 0.001, no origin point in the AUPRC). Tables
scored before 2026-10-15 need ``--min-tau 0 --anchor-origin``.

Regeneration command (rung outputs are the byte-identical re-runs of the §4
factorial — one pipeline run per sub-directory of ``--run-dir``)::

    uv run python validation/ablation_pr_decomposition.py \\
        --t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz \\
        --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \\
        --run-dir <root with one sub-directory per factorial rung>
"""

from __future__ import annotations

import importlib.util
import inspect
import sys
from collections import defaultdict
from pathlib import Path

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


MIN_TAU = getattr(rs, "CAFA_MIN_TAU", 0.0)


def _metric_kwargs(min_tau: float, anchor_origin: bool) -> dict:
    """Convention kwargs ``resampling.panel_metrics`` accepts (none, if older)."""
    params = inspect.signature(rs.panel_metrics).parameters
    wanted = {"min_tau": min_tau, "anchor_origin": anchor_origin}
    return {k: v for k, v in wanted.items() if k in params}


def f_max_index(panel, min_tau: float = MIN_TAU) -> int | None:
    """Index of the F_max cutoff among the thresholds ``>= min_tau``.

    The same rule ``resampling.panel_metrics`` applies; ``None`` when no cutoff
    is reachable or every reachable cutoff has F = 0 (``panel_metrics`` then
    reports no F_max cutoff either, so there is no operating point to read).
    """
    reachable = np.nonzero(panel.thresholds >= min_tau)[0]
    if not len(reachable):
        return None
    f = rs.panel_curve(panel).f
    best = int(reachable[np.argmax(f[reachable])])
    return best if f[best] > 0 else None


def pr_at_fmax(panel, min_tau: float = MIN_TAU) -> dict[str, float]:
    """Precision/recall/coverage at the F_max operating point of a panel.

    Reads the full-cohort CAFA curves off ``resampling.panel_curve`` — the
    exact vectors ``panel_metrics`` maximises — takes ``argmax F`` over the
    cutoffs ``>= min_tau`` and reports PPV and sensitivity at that index.
    """
    empty = {
        "f_max": 0.0,
        "f_max_tau": 0.0,
        "ppv_at_fmax": 0.0,
        "sensitivity_at_fmax": 0.0,
        "coverage_at_fmax": 0.0,
    }
    if panel.n_proteins == 0:
        return empty
    best = f_max_index(panel, min_tau)
    if best is None:
        return empty
    curve = rs.panel_curve(panel)
    f_max = float(curve.f[best])
    return {
        "f_max": f_max,
        "f_max_tau": float(panel.thresholds[best]) if f_max > 0 else 0.0,
        "ppv_at_fmax": float(curve.precision[best]),
        "sensitivity_at_fmax": float(curve.recall[best]),
        "coverage_at_fmax": float(curve.coverage[best]),
    }


def _interp_upper(x: np.ndarray, y: np.ndarray, at: float) -> float:
    """``y`` at ``x = at`` on the upper envelope of the points, NaN outside."""
    best: dict[float, float] = {}
    for xi, yi in zip(x.tolist(), y.tolist()):
        if xi not in best or yi > best[xi]:
            best[xi] = yi
    xs = np.array(sorted(best))
    if xs.size == 0 or at < xs[0] or at > xs[-1]:
        return float("nan")
    return float(np.interp(at, xs, np.array([best[v] for v in xs])))


def reference_operating_point(panel, min_tau: float = MIN_TAU) -> tuple[float, float]:
    """``(sensitivity, n_predictions)`` of a panel at its F_max threshold."""
    best = f_max_index(panel, min_tau)
    if best is None:
        return float("nan"), float("nan")
    curve = rs.panel_curve(panel)
    return float(curve.recall[best]), float(panel.n_pred[:, best].sum())


def ppv_at_sensitivity(panel, sensitivity: float) -> float:
    """PPV of a panel's own curve at ``sensitivity`` (upper envelope, interpolated)."""
    curve = rs.panel_curve(panel)
    return _interp_upper(curve.recall, curve.precision, sensitivity)


def f_at_budget(panel, n_predictions: float) -> float:
    """F of a panel's own curve at ``n_predictions`` total predictions."""
    curve = rs.panel_curve(panel)
    totals = panel.n_pred.sum(axis=0)
    return _interp_upper(totals, curve.f, n_predictions)


def main() -> int:  # pragma: no cover - I/O wiring
    import argparse

    import pandas as pd
    from loguru import logger

    logger.remove()
    logger.add(sys.stderr, level="INFO")

    if str(_ROOT) not in sys.path:
        sys.path.insert(0, str(_ROOT))

    from src.domain_annotation_parser import DomainAnnotationParser
    from src.goa_parser import EXPERIMENTAL_EVIDENCE, GOAParser, parse_goa
    from src.ontology_processor import OntologyProcessor

    parser = argparse.ArgumentParser(
        description="Sensitivity/PPV decomposition at the F_max operating point "
        "for every §4 ablation rung."
    )
    parser.add_argument(
        "--t0-gaf", type=Path, default=Path("data/raw/goa_archive/goa_human.gaf.205.gz")
    )
    parser.add_argument(
        "--t1-gaf", type=Path, default=Path("data/raw/goa_annotations/goa_human.gaf.gz")
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=Path("results/ablation"),
        help="Root holding one sub-directory per current-code factorial run",
    )
    parser.add_argument(
        "--interpro", type=Path, default=Path("data/interim/protein2ipr_human.dat.gz")
    )
    parser.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    parser.add_argument(
        "--committed-metrics",
        type=Path,
        default=_HERE / "ablation_metrics.tsv",
        help="The committed §4 metrics table the recomputation must reproduce",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=_HERE / "ablation_pr_decomposition.tsv",
    )
    parser.add_argument(
        "--reference-rung",
        default="supra_input_relative",
        help="Rung whose F_max operating point fixes the sensitivity and the "
        "prediction budget the other rungs are read at "
        "(default: supra_input_relative)",
    )
    parser.add_argument(
        "--min-tau",
        type=float,
        default=MIN_TAU,
        help="Lowest cutoff F_max is read at (default: %(default)s, the CAFA "
        "convention). Tables scored before 2026-10-15 used the whole sweep: 0",
    )
    parser.add_argument(
        "--anchor-origin",
        action="store_true",
        help="Keep the zero-coverage points in the AUPRC trapezoid (the "
        "pre-2026-10-15 convention; default: leave them out, as CAFA does)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-6,
        help="Sanity gate: max |recomputed - committed| f_max/auprc allowed",
    )
    args = parser.parse_args()
    ic_floors = [0.0, 2.0, 4.0]
    if args.reference_rung not in {r.name for r in abl.LADDER}:
        parser.error(f"--reference-rung {args.reference_rung} is not a ladder rung")

    for path in (
        args.t0_gaf,
        args.t1_gaf,
        args.interpro,
        args.go_ontology,
        args.committed_metrics,
    ):
        if not path.exists():
            logger.error(f"Missing required input: {path}")
            return 1

    # -------------------------------------------------- assembly (as ablation) --
    logger.info("Loading GO ontology...")
    processor = OntologyProcessor(args.go_ontology)
    get_ancestors = processor.get_ancestors
    term_aspect = tb.build_term_aspect(processor)

    logger.info("Parsing t0 (training) GOA — all non-IEA evidence...")
    t0_map = parse_goa(args.t0_gaf, evidence_filter="manual")

    logger.info("Parsing t1 (test) GOA — experimental evidence only...")
    t1_exp_map = GOAParser(
        evidence_codes=EXPERIMENTAL_EVIDENCE, aspects={"P", "F", "C"}
    ).parse_gaf_file(args.t1_gaf)
    t0_map = tb.remap_retired_ids(
        t0_map, processor.alt_id_map, processor.replaced_by_map
    )
    t1_exp_map = tb.remap_retired_ids(
        t1_exp_map, processor.alt_id_map, processor.replaced_by_map
    )

    logger.info("Parsing domain architectures...")
    dom_parser = DomainAnnotationParser(max_supra_domain_length=3, min_domain_length=10)
    architectures = dom_parser.parse_protein2ipr_file(args.interpro)
    protein_domains: dict[str, list[str]] = {}
    for protein, arch in architectures.items():
        domains = list(arch.single_domains) + list(arch.supra_domains)
        if domains:
            protein_domains[protein] = domains

    logger.info("Building the CAFA no-knowledge benchmark...")
    benchmark = tb.build_nk_benchmark_by_aspect(
        t0_map,
        t1_exp_map,
        term_aspect,
        get_ancestors,
        predictable_proteins=set(protein_domains),
    )
    eval_proteins = {p for aspect in benchmark.values() for p in aspect}
    for aspect, truths in benchmark.items():
        logger.info(f"  {aspect}: {len(truths):,} no-knowledge benchmark proteins")
    if not eval_proteins:
        logger.error("Empty benchmark.")
        return 1

    logger.info("Computing information content from t0...")
    ic = tb.information_content(t0_map, get_ancestors)
    n_t0 = len(t0_map)
    freq_counts: dict[str, int] = defaultdict(int)
    for terms in t0_map.values():
        for t in tb.propagate_terms(terms, get_ancestors):
            freq_counts[t] += 1
    term_freq = {t: c / n_t0 for t, c in freq_counts.items()} if n_t0 else {}

    # ----------------------------------------------------- transfer per rung --
    eval_domains = {p: protein_domains[p] for p in eval_proteins}
    methods: dict[str, dict] = {}
    for rung in abl.LADDER:
        path = abl.rung_prediction_file(args.run_dir, rung)
        if not path.exists():
            logger.error(f"[{rung.name}] missing predictions: {path}")
            return 1
        scores_q = tb.remap_retired_scores(
            tb.load_domain_go_scores(
                path,
                "q_value" if rung.kind == "propagated" else "adj_p_value",
                neg_log10=True,
            ),
            processor.alt_id_map,
            processor.replaced_by_map,
            term_aspect,
        )
        methods[rung.name] = tb.transfer_predictions_pscore(
            eval_domains, scores_q, get_ancestors
        )
        logger.info(f"  [{rung.name}] {len(scores_q):,} domains scored")
    methods["naive"] = tb.naive_predictions(eval_proteins, term_freq)

    aspect_preds = {
        name: {
            aspect: tb.restrict_to_aspect(preds, aspect, term_aspect)
            for aspect in ("BP", "MF", "CC")
        }
        for name, preds in methods.items()
    }

    # ------------------------------------------------- evaluate + decompose --
    committed = pd.read_csv(args.committed_metrics, sep="\t")
    committed_lookup = {
        (r["aspect"], float(r["min_ic"]), r["method"]): r
        for _, r in committed.iterrows()
    }

    rows: list[dict] = []
    mismatches: list[str] = []
    max_fmax_delta = 0.0
    max_auprc_delta = 0.0
    for aspect in ("BP", "MF", "CC"):
        if not benchmark[aspect]:
            continue
        for min_ic in ic_floors:
            true_a = tb.filter_by_ic(benchmark[aspect], ic, min_ic)
            if not true_a:
                continue
            panels = {}
            for name in methods:
                pred_f = tb.filter_by_ic(aspect_preds[name][aspect], ic, min_ic)
                taus = tb._candidate_thresholds(pred_f)
                panels[name] = rs.build_panel(pred_f, true_a, ic, taus)
            ref_sensitivity, ref_budget = reference_operating_point(
                panels[args.reference_rung], args.min_tau
            )
            for name, panel in panels.items():
                observed = rs.panel_metrics(
                    panel, **_metric_kwargs(args.min_tau, args.anchor_origin)
                )
                decomp = pr_at_fmax(panel, args.min_tau)
                if abs(decomp["f_max"] - observed["f_max"]) > 1e-12:
                    logger.error(
                        f"[{aspect} IC>={min_ic:g}] {name}: pr_at_fmax disagrees "
                        f"with panel_metrics ({decomp['f_max']} vs "
                        f"{observed['f_max']}) — internal bug, aborting."
                    )
                    return 1

                key = (aspect, min_ic, name)
                com = committed_lookup.get(key)
                if com is None:
                    mismatches.append(f"{key}: missing from committed metrics")
                    continue
                d_f = abs(observed["f_max"] - com["f_max"])
                d_a = abs(observed["auprc"] - com["auprc"])
                max_fmax_delta = max(max_fmax_delta, d_f)
                max_auprc_delta = max(max_auprc_delta, d_a)
                if d_f > args.tolerance or d_a > args.tolerance:
                    mismatches.append(
                        f"{key}: f_max {observed['f_max']:.12f} vs committed "
                        f"{com['f_max']:.12f} (Δ={d_f:.3e}); auprc "
                        f"{observed['auprc']:.12f} vs {com['auprc']:.12f} "
                        f"(Δ={d_a:.3e})"
                    )
                if observed["n_eval_proteins"] != int(com["n_eval_proteins"]):
                    mismatches.append(
                        f"{key}: n_eval_proteins {observed['n_eval_proteins']} vs "
                        f"committed {int(com['n_eval_proteins'])}"
                    )

                rows.append(
                    {
                        "aspect": aspect,
                        "min_ic": min_ic,
                        "rung": name,
                        "n_eval_proteins": observed["n_eval_proteins"],
                        "f_max": observed["f_max"],
                        "f_max_tau": observed["f_max_tau"],
                        "ppv_at_fmax": decomp["ppv_at_fmax"],
                        "sensitivity_at_fmax": decomp["sensitivity_at_fmax"],
                        "coverage_at_fmax": decomp["coverage_at_fmax"],
                        "reference_rung": args.reference_rung,
                        "ref_sensitivity": ref_sensitivity,
                        "ref_n_predictions": ref_budget,
                        "ppv_at_ref_sensitivity": ppv_at_sensitivity(
                            panel, ref_sensitivity
                        ),
                        "f_at_budget": f_at_budget(panel, ref_budget),
                        "committed_f_max": com["f_max"],
                        "f_max_delta": observed["f_max"] - com["f_max"],
                    }
                )
                logger.info(
                    f"  [{aspect} IC>={min_ic:g}] {name:22s} "
                    f"F_max={observed['f_max']:.4f} "
                    f"PPV={decomp['ppv_at_fmax']:.4f} "
                    f"sens={decomp['sensitivity_at_fmax']:.4f} "
                    f"cov={decomp['coverage_at_fmax']:.2f} (Δf={d_f:.2e})"
                )

    # ------------------------------------------------------------ sanity gate --
    logger.info(
        f"Sanity gate: max |Δf_max|={max_fmax_delta:.3e}, "
        f"max |Δauprc|={max_auprc_delta:.3e} vs {args.committed_metrics}"
    )
    if mismatches:
        logger.error(
            f"SANITY GATE FAILED — {len(mismatches)} cell(s) disagree with the "
            f"committed metrics beyond {args.tolerance:g}. Not writing output."
        )
        for line in mismatches:
            logger.error(f"  {line}")
        return 1

    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    logger.info(f"✓ PR decomposition ({len(rows)} rows): {args.output}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
