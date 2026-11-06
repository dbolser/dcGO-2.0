"""Reliability of the domain-unit scores by the support of the winning unit.

A protein's domain-unit score for a term is the largest θ(u, t) over its
units. When a unit is carried by only one to three annotated proteins, its θ
is a noisy estimate, and taking the maximum over many units picks the most
over-estimated ones (a winner's curse). This post-hoc analysis rebuilds the
frozen domain-unit arm on a split (same inputs and flags as the frozen run),
records for every prediction the number of annotated carriers of the unit
that won the maximum, and tabulates precision at IC >= 2 by score band and
support band, with the number of distinct proteins and a 95% interval from
1,000 resamples of proteins (all of a protein's predictions move together).
Writes ``<run>/reliability_by_support.tsv``.

    uv run --group projection python paper/reliability_by_support.py dev
    uv run --group projection python paper/reliability_by_support.py final
    uv run --group projection python paper/reliability_by_support.py final humanonly

Every run also writes ``<run>/reliability_contrasts.tsv``: paired
protein-bootstrap differences in precision (units with 1–3 carriers minus
units with more than 50; with the background, other-species-only units
minus units with a human carrier) within each score band at 0.5 and above,
and per support band the observed precision minus the mean score
(a calibration gap). With the all-species background it also writes
``<run>/reliability_by_origin.tsv``: the same precision by score band and by
where the winning unit's carriers come from (human only, human and other
species, other species only). ``humanonly`` trains on human annotations
alone (as ``final_humanonly``) and writes into that run's directory.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
_spec = importlib.util.spec_from_file_location(
    "projection_benchmark", ROOT / "validation" / "projection_benchmark.py"
)
pb = importlib.util.module_from_spec(_spec)
sys.modules["projection_benchmark"] = pb
_spec.loader.exec_module(pb)

from paper.projection_training_sizes import SPLITS  # noqa: E402
from src import knowledge_projection as kp  # noqa: E402
from validation.knowledge_yield import bootstrap_weights  # noqa: E402

OUT_DIR = {
    ("dev", False): "dev2018_frozen",
    ("final", False): "final",
    ("dev", True): "dev2018_humanonly",
    ("final", True): "final_humanonly",
}
ORIGINS = ["human only", "human and other species", "other species only"]


def _weighted(g: pd.DataFrame, weights: np.ndarray, column: str) -> np.ndarray:
    """Per-resample mean of ``column`` over a group, proteins resampled."""
    per = g.groupby("protein")[column].agg(["sum", "size"])
    w = weights[:, per.index.to_numpy()]
    total = w @ per["size"].to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        return (w @ per["sum"].to_numpy()) / total


def _interval(x: np.ndarray) -> tuple[float, float]:
    x = x[np.isfinite(x)]
    return (float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5)))


def contrast(a: pd.DataFrame, b: pd.DataFrame, weights: np.ndarray) -> dict:
    """Precision of group ``a`` minus group ``b``, paired over the same resamples."""
    d = _weighted(a, weights, "correct") - _weighted(b, weights, "correct")
    lo, hi = _interval(d)
    return {
        "n_a": len(a),
        "proteins_a": a.protein.nunique(),
        "precision_a": float(a.correct.mean()),
        "n_b": len(b),
        "proteins_b": b.protein.nunique(),
        "precision_b": float(b.correct.mean()),
        "difference": float(a.correct.mean() - b.correct.mean()),
        "ci_lo": lo,
        "ci_hi": hi,
    }


def calibration_gap(g: pd.DataFrame, weights: np.ndarray) -> dict:
    """Observed precision minus mean score, with a protein-bootstrap interval."""
    d = _weighted(g, weights, "correct") - _weighted(g, weights, "score")
    lo, hi = _interval(d)
    return {
        "n": len(g),
        "proteins": g.protein.nunique(),
        "precision": float(g.correct.mean()),
        "mean_score": float(g.score.mean()),
        "gap": float(g.correct.mean() - g.score.mean()),
        "ci_lo": lo,
        "ci_hi": hi,
    }


def cell_stats(g: pd.DataFrame, weights: np.ndarray) -> dict:
    """Precision of one cell with a protein-bootstrap interval."""
    per_protein = g.groupby("protein").correct.agg(["sum", "size"])
    w = weights[:, per_protein.index.to_numpy()]
    hits = w @ per_protein["sum"].to_numpy()
    total = w @ per_protein["size"].to_numpy()
    reps = hits[total > 0] / total[total > 0]
    return {
        "predictions": len(g),
        "proteins": len(per_protein),
        "correct": int(g.correct.sum()),
        "precision": float(g.correct.mean()),
        "ci_lo": float(np.percentile(reps, 2.5)),
        "ci_hi": float(np.percentile(reps, 97.5)),
    }


N_BOOT = 1000
SUPPORT_BINS = [0, 3, 10, 50, np.inf]
SUPPORT_LABELS = ["1–3", "4–10", "11–50", "> 50"]
SCORE_BINS = [0.0, 0.3, 0.5, 0.7, 0.9, 1.0001]
SCORE_LABELS = ["< 0.3", "0.3–0.5", "0.5–0.7", "0.7–0.9", "≥ 0.9"]


def main() -> int:
    split = sys.argv[1]
    human_only = len(sys.argv) > 2 and sys.argv[2] == "humanonly"
    args = pb.parse_args(
        SPLITS[split] + ["--extra-evidence-filter", "experimental", "--strength", "0.5"]
    )
    inputs = pb.ab.load_benchmark_inputs(
        args.t0_gaf,
        args.t1_gaf,
        args.interpro,
        args.go_ontology,
        "interpro",
        args.evidence_filter,
        cohort_interpro=args.cohort_interpro,
    )
    families = pb.unit_families(
        inputs.protein_domains, args.interpro, set(inputs.protein_domains)
    )
    human = set(inputs.protein_domains)
    if human_only:
        merged = {p: set(ts) for p, ts in inputs.t0_map.items()}
    else:
        merged = pb.add_background(args, inputs, families)
    trainings = kp.build_aspect_trainings(
        merged,
        inputs.get_ancestors,
        inputs.term_aspect,
        pb.ASPECTS,
        exclude={a: list(inputs.benchmark[a]) for a in pb.ASPECTS},
    )
    rows, origin_rows, contrast_rows = [], [], []
    for aspect in pb.ASPECTS:
        truth = pb.tb.filter_by_ic(inputs.benchmark[aspect], inputs.ic, 2.0)
        queries = sorted(truth)
        training = trainings[aspect]
        vocab = kp.UnitVocabulary(pb.combine_families(families, ["ipr"]))
        train_inc = vocab.incidence(training.proteins)
        evidence = kp.unit_term_evidence(train_inc, training)
        is_human = np.array([p in human for p in training.proteins])
        human_carriers = np.asarray(train_inc[is_human].sum(axis=0)).ravel()
        theta = evidence.posterior(args.strength[0])
        scores, arg = kp.max_product(
            vocab.incidence(queries), theta, return_argmax=True
        )
        s, a = scores.tocoo(), arg.tocsr()
        unit = np.asarray(a[s.row, s.col]).ravel().astype(int) - 1
        terms = np.array(training.terms)[s.col]
        keep = np.array([inputs.ic.get(t, 0.0) >= 2.0 for t in terms])
        hit = np.array(
            [t in truth[queries[r]] for r, t in zip(s.row[keep], terms[keep])]
        )
        frame = pd.DataFrame(
            {
                "protein": s.row[keep],
                "score": s.data[keep],
                "support": evidence.n_carriers[unit[keep]],
                "human_support": human_carriers[unit[keep]],
                "correct": hit,
            }
        )
        frame["support_band"] = pd.cut(
            frame.support, SUPPORT_BINS, labels=SUPPORT_LABELS
        )
        frame["score_band"] = pd.cut(
            frame.score, SCORE_BINS, labels=SCORE_LABELS, right=False
        )
        frame["origin"] = np.select(
            [frame.human_support == frame.support, frame.human_support > 0],
            ORIGINS[:2],
            ORIGINS[2],
        )
        weights = bootstrap_weights(len(queries), N_BOOT, seed=0)
        for keys, g in frame.groupby(["score_band", "support_band"], observed=True):
            rows.append(
                {"aspect": aspect, "score_band": keys[0], "support_band": keys[1]}
                | cell_stats(g, weights)
            )
        for keys, g in frame.groupby(["score_band", "origin"], observed=True):
            origin_rows.append(
                {"aspect": aspect, "score_band": keys[0], "origin": keys[1]}
                | cell_stats(g, weights)
            )
        top = frame[frame.score >= 0.5]
        for band in ["0.5–0.7", "0.7–0.9", "≥ 0.9", "≥ 0.5"]:
            g = top if band == "≥ 0.5" else top[top.score_band == band]
            pairs = {
                "support 1–3 − support > 50": (
                    g[g.support_band == "1–3"],
                    g[g.support_band == "> 50"],
                ),
            }
            if not human_only:
                pairs["other species only − any human carrier"] = (
                    g[g.origin == ORIGINS[2]],
                    g[g.origin != ORIGINS[2]],
                )
            for name, (a, b) in pairs.items():
                if a.protein.nunique() >= 2 and b.protein.nunique() >= 2:
                    contrast_rows.append(
                        {"aspect": aspect, "kind": "contrast", "score_band": band}
                        | {"comparison": name}
                        | contrast(a, b, weights)
                    )
        for band in ["1–3", "4–10", "11–50", "> 50"]:
            g = top[top.support_band == band]
            if g.protein.nunique() >= 2:
                contrast_rows.append(
                    {"aspect": aspect, "kind": "calibration", "score_band": "≥ 0.5"}
                    | {"comparison": f"support {band}: precision − mean score"}
                    | calibration_gap(g, weights)
                )
        print(aspect, len(frame), "predictions at IC >= 2")
    out_dir = ROOT / "validation/projection" / OUT_DIR[split, human_only]
    pd.DataFrame(rows).to_csv(
        out_dir / "reliability_by_support.tsv", sep="\t", index=False
    )
    pd.DataFrame(contrast_rows).to_csv(
        out_dir / "reliability_contrasts.tsv", sep="\t", index=False
    )
    print(pd.DataFrame(contrast_rows).round(3).to_string(index=False))
    if not human_only:
        pd.DataFrame(origin_rows).to_csv(
            out_dir / "reliability_by_origin.tsv", sep="\t", index=False
        )
        print(pd.DataFrame(origin_rows).to_string(index=False))
    print(pd.DataFrame(rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
