"""Does the all-species background hurt domain units through other-species-only units?

Post-hoc check of the explanation in Section 5.3. On the test split, builds
the frozen domain-unit arm three ways:

- ``human only``: trained on human t0 annotations alone (as ``final_humanonly``);
- ``+ background``: the frozen arm (human + all-species experimental);
- ``+ background, no other-species-only units``: the frozen arm with θ set to
  0 for every unit none of whose annotated carriers is a human protein.

Reports KY50 and pooled AUPRC at IC >= 2 with paired protein-bootstrap
intervals (1,000 resamples, the same resamples for every arm) for each arm
and for its differences from the frozen arm and from the human-only arm.
Writes ``validation/projection/final/background_origin_ablation.tsv`` and
the human-only arm's precision–yield curve
(``background_origin_curves.tsv``, for Figure 1).

    uv run --group projection python paper/background_origin_ablation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from paper.projection_training_sizes import SPLITS, pb  # noqa: E402
from src import knowledge_projection as kp  # noqa: E402
from validation import knowledge_yield as ky  # noqa: E402

N_BOOT = 1000
FROZEN = "+ background"
HUMAN = "human only"


def domain_scores(training, families, queries, strength, keep_unit=None):
    """``{protein: {term: θ-max score}}`` for the domain-unit arm."""
    vocab = kp.UnitVocabulary(pb.combine_families(families, ["ipr"]))
    train_inc = vocab.incidence(training.proteins)
    theta = kp.unit_term_evidence(train_inc, training).posterior(strength)
    if keep_unit is not None:
        theta = sp.diags(keep_unit(train_inc).astype(float)) @ theta
    scores = kp.max_product(vocab.incidence(queries), sp.csr_matrix(theta)).tocoo()
    out: dict[str, dict[str, float]] = {}
    for r, c, v in zip(scores.row, scores.col, scores.data):
        out.setdefault(queries[r], {})[training.terms[c]] = float(v)
    return out


def main() -> int:
    args = pb.parse_args(
        SPLITS["final"]
        + ["--extra-evidence-filter", "experimental", "--strength", "0.5"]
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
    human = set(inputs.protein_domains)
    exclude = {a: list(inputs.benchmark[a]) for a in pb.ASPECTS}
    human_families = pb.unit_families(inputs.protein_domains, args.interpro, human)
    human_train = kp.build_aspect_trainings(
        {p: set(t) for p, t in inputs.t0_map.items()},
        inputs.get_ancestors,
        inputs.term_aspect,
        pb.ASPECTS,
        exclude=exclude,
    )
    bg_families = pb.unit_families(inputs.protein_domains, args.interpro, human)
    merged = pb.add_background(args, inputs, bg_families)
    bg_train = kp.build_aspect_trainings(
        merged, inputs.get_ancestors, inputs.term_aspect, pb.ASPECTS, exclude=exclude
    )
    rows, curves = [], []
    for aspect in pb.ASPECTS:
        truth = pb.tb.filter_by_ic(inputs.benchmark[aspect], inputs.ic, 2.0)
        order = sorted(truth)
        sizes = np.array([len(truth[p]) for p in order], dtype=float)
        bg = bg_train[aspect]
        is_human = np.array([p in human for p in bg.proteins])

        def has_human_carrier(inc, is_human=is_human):
            return np.asarray(inc[is_human].sum(axis=0)).ravel() > 0

        arms = {
            HUMAN: domain_scores(
                human_train[aspect], human_families, order, args.strength[0]
            ),
            FROZEN: domain_scores(bg, bg_families, order, args.strength[0]),
            "+ background, no other-species-only units": domain_scores(
                bg, bg_families, order, args.strength[0], has_human_carrier
            ),
        }
        rankings = {
            name: ky.pool(
                truth,
                pb.tb.filter_by_ic(preds, inputs.ic, 2.0),
                inputs.ic,
                proteins=order,
            )
            for name, preds in arms.items()
        }
        reps = ky.bootstrap(rankings, sizes, n_replicates=N_BOOT, seed=0)
        for name, ranking in rankings.items():
            point = ky.summarise(ranking)
            row = {"aspect": aspect, "arm": name}
            for metric in ("ky50", "pooled_auprc"):
                rep = reps[name][metric]
                row |= {
                    metric: point[metric],
                    f"{metric}_ci_lo": float(np.percentile(rep, 2.5)),
                    f"{metric}_ci_hi": float(np.percentile(rep, 97.5)),
                }
                for ref, tag in ((FROZEN, "frozen"), (HUMAN, "human")):
                    diff = rep - reps[ref][metric]
                    row |= {
                        f"{metric}_diff_vs_{tag}": point[metric]
                        - ky.summarise(rankings[ref])[metric],
                        f"{metric}_diff_vs_{tag}_ci_lo": float(
                            np.percentile(diff, 2.5)
                        ),
                        f"{metric}_diff_vs_{tag}_ci_hi": float(
                            np.percentile(diff, 97.5)
                        ),
                    }
            rows.append(row)
        human_rank = rankings[HUMAN]
        tp, pp, _ = ky._curves(human_rank, None)
        thr = human_rank.scores[human_rank.group_end - 1]
        prec = np.where(pp > 0, tp / np.maximum(pp, 1e-300), 0.0)
        keep = prec >= 0.1  # the plotted range; the long tail is dropped
        for t, y, x, q in zip(thr[keep], tp[keep], pp[keep], prec[keep]):
            curves.append(
                {
                    "aspect": aspect,
                    "min_ic": 2.0,
                    "method": "domain units, human training only",
                    "threshold": t,
                    "tp": y,
                    "predicted": x,
                    "precision": q,
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(
        ROOT / "validation/projection/final/background_origin_ablation.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(curves).to_csv(
        ROOT / "validation/projection/final/background_origin_curves.tsv",
        sep="\t",
        index=False,
    )
    print(out.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
