#!/usr/bin/env python3
"""Domain-centric evaluation of dcGO associations against InterPro2GO.

Companion to the protein-centric temporal benchmark (`temporal_benchmark.py`).
Where that scores *proteins*, this scores the **(domain, GO) associations
themselves** — the unit dcGO actually produces — against the curated InterPro2GO
map, with GO-DAG propagation and restricted to the domains present in both sets
(the §1 reference frame). It is the right instrument for judging the *relative
(parental-background) inference*, whose job is to raise the specificity/precision
of the association set — an effect the protein-centric F_max is largely blind to.

Pass one or more labelled association TSVs and they are compared side by side:

    uv run python validation/domain_centric_eval.py \
        --predictions base=results_t0_2021/domain_go_associations_significant.tsv \
        --predictions relative=results_t0_2021/domain_go_associations_relative.tsv

Metrics (on the shared, propagated domain space):
  * precision   = recovered / predicted   (fraction of our pairs that are curated;
                  a *lower bound* — InterPro2GO is incomplete, so genuinely-novel
                  pairs count against it. Useful for *comparing* configurations.)
  * recall      = recovered / reference    (coverage of InterPro2GO)
  * F1          = harmonic mean

Note: InterPro2GO covers single InterPro entries only, so supra-domains are
dropped by the shared-domain restriction — this evaluates the single-domain
associations. The reference here is the *current* InterPro2GO; swapping in a
dated (e.g. 2021) InterPro2GO would make it a temporal domain-centric test.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import pandas as pd
from loguru import logger

logger.remove()
logger.add(sys.stderr, level="INFO")

# Reuse the §1 helpers (validation/ is not a package; load by path).
_VR = Path(__file__).resolve().parent / "validate_results.py"
_spec = importlib.util.spec_from_file_location("validate_results", _VR)
vr = importlib.util.module_from_spec(_spec)
sys.modules["validate_results"] = vr
_spec.loader.exec_module(vr)


def evaluate_one(
    predictions_file: Path,
    reference: set,
    get_ancestors,
) -> dict:
    """Domain-centric precision/recall of one association set vs InterPro2GO."""
    df = pd.read_csv(predictions_file, sep="\t")
    pred_pairs = set(zip(df["domain"], df["go_term"]))

    # Shared single-domain space, then propagate both sides to ancestor closure.
    shared = {d for d, _ in pred_pairs} & {d for d, _ in reference}
    pred = vr.propagate_pairs({p for p in pred_pairs if p[0] in shared}, get_ancestors)
    ref = vr.propagate_pairs({p for p in reference if p[0] in shared}, get_ancestors)
    m = vr.compute_metrics(pred, ref, str(predictions_file))
    return {
        "n_associations": len(df),
        "shared_domains": len(shared),
        "predicted_pairs": m["n_predictions"],
        "recovered": m["recovered"],
        "precision_lb": m["precision_lower_bound"],
        "recall": m["reference_coverage"],
        "f1": (
            2
            * m["precision_lower_bound"]
            * m["reference_coverage"]
            / (m["precision_lower_bound"] + m["reference_coverage"])
            if (m["precision_lower_bound"] + m["reference_coverage"]) > 0
            else 0.0
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--predictions",
        action="append",
        required=True,
        metavar="LABEL=PATH",
        help="Labelled association TSV (repeatable), e.g. base=path.tsv",
    )
    ap.add_argument(
        "--reference", type=Path, default=Path("data/raw/interpro2go/interpro2go")
    )
    ap.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    ap.add_argument("--output-dir", type=Path, default=Path("validation"))
    args = ap.parse_args()

    if not args.reference.exists():
        alt = Path("data/raw/interpro_mappings/interpro2go")
        if alt.exists():
            args.reference = alt
    for label_path in args.predictions:
        if "=" not in label_path:
            logger.error(f"--predictions must be LABEL=PATH, got: {label_path}")
            return 1
    if not args.reference.exists() or not args.go_ontology.exists():
        logger.error(
            f"Missing reference ({args.reference}) or ontology ({args.go_ontology}). "
            "Download with: uv run python scripts/download_data.py "
            "--datasets interpro2go --datasets go_ontology"
        )
        return 1

    reference = vr.parse_interpro2go(args.reference)
    get_ancestors = vr.build_get_ancestors(args.go_ontology)

    rows = []
    for label_path in args.predictions:
        label, path = label_path.split("=", 1)
        p = Path(path)
        if not p.exists():
            logger.error(f"Predictions file not found: {p}")
            return 1
        logger.info(f"Evaluating '{label}': {p}")
        res = evaluate_one(p, reference, get_ancestors)
        res["config"] = label
        rows.append(res)
        logger.info(
            f"  precision(lb)={res['precision_lb']:.3f}  recall={res['recall']:.3f}  "
            f"F1={res['f1']:.3f}  (predicted pairs={res['predicted_pairs']:,})"
        )

    df = pd.DataFrame(rows)[
        [
            "config",
            "n_associations",
            "shared_domains",
            "predicted_pairs",
            "recovered",
            "precision_lb",
            "recall",
            "f1",
        ]
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir / "domain_centric_metrics.tsv"
    df.to_csv(out, sep="\t", index=False)
    logger.info(f"✓ Saved: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
