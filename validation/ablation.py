#!/usr/bin/env python3
"""Component ablation + permutation null + bootstrap CIs — VALIDATION_PLAN §4.

Closes three P0 publication blockers from `ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md`
in one driver, because they all need the same expensive scaffolding (the §2
temporal split, its no-knowledge cohort, and its propagated truth):

**(a) The ablation ladder.** Which parts of the method earn their place?

    ==================  =====================================================
    rung                what it adds
    ==================  =====================================================
    ``single``          single InterPro domains only
    ``supra``           + supra-domains (contiguous combinations, len <= 3)
    ``supra_input``     + True Path propagation of input annotations
    ``supra_relative``  + relative inference
    ``supra_output``    + True Path propagation of inferred associations
    ``full``            + all three hierarchy stages
    ==================  =====================================================

**(b) A permutation null**, not one shuffle: ``--n-permutations`` seeded
re-labellings of the domain -> term map, summarised as a null distribution with a
percentile interval and an empirical p-value.

**(c) Protein-level bootstrap CIs**, including the **paired** difference between
rungs and against the naive baseline — resampling the benchmark proteins once per
replicate and recomputing every method on that resample.

Plus the two P1 reporting items that fall out of it: the full selection-stage
count ledger (including how the cohort changes with the IC floor), and prediction
coverage next to F_max.

How the rungs are produced
--------------------------
Fisher tests are independent per (feature, term) cell, so the rungs share work
without sharing statistics:

* ``single`` is its own pipeline run (``--disable-supra-domains``) because its BH
  hypothesis family is genuinely smaller — 3.1e8 tests, not 1.6e9 — and the FDR
  cut therefore differs.
* Every hierarchy-stage combination is a complete pipeline run. This is
  essential because input propagation changes the test universe and relative
  inference is combined before BH correction.
* Output propagation rungs are read from the pipeline's
  ``domain_go_annotations_propagated.tsv``; all other rungs use the significant
  association table.

Honest caveat, stated once here and again in ``VALIDATION_PLAN.md``: the True
Path rungs are scored on ``q_value`` because the propagated output carries no
``p_value`` column. To keep the ladder comparable, **every** rung's primary score
is ``-log10(q)``; the ``-log10(p)`` variant of the non-TPR rungs is reported
alongside as a sensitivity check (``score`` column of the metrics TSV).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

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


tb = _load_sibling("temporal_benchmark")
rs = _load_sibling("resampling")


# --------------------------------------------------------------------------- #
# The ladder                                                                   #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Rung:
    """One step of the ablation ladder."""

    name: str
    label: str
    adds: str
    #: sub-directory of --run-dir holding this rung's pipeline output
    run_dir: str
    #: which file in it: the direct associations or the True-Path-propagated ones
    kind: str  # "associations" | "propagated"
    #: rung this one is built on (for the paired "did this component help?" test)
    parent: str | None


LADDER: tuple[Rung, ...] = (
    Rung("single", "single domains", "-", "single", "associations", None),
    Rung(
        "supra", "+ supra-domains", "supra-domains", "supra", "associations", "single"
    ),
    Rung(
        "supra_input",
        "+ input propagation",
        "input propagation",
        "supra_input",
        "associations",
        "supra",
    ),
    Rung(
        "supra_relative",
        "+ relative inference",
        "relative inference",
        "supra_relative",
        "associations",
        "supra",
    ),
    Rung(
        "supra_output",
        "+ output propagation",
        "output propagation",
        "supra_output",
        "propagated",
        "supra",
    ),
    Rung(
        "supra_input_relative",
        "input + relative",
        "input propagation + relative inference",
        "supra_input_relative",
        "associations",
        "supra_input",
    ),
    Rung(
        "supra_input_output",
        "input + output",
        "input + output propagation",
        "supra_input_output",
        "propagated",
        "supra_input",
    ),
    Rung(
        "supra_relative_output",
        "relative + output",
        "relative inference + output propagation",
        "supra_relative_output",
        "propagated",
        "supra_relative",
    ),
    Rung(
        "full",
        "full current method",
        "input propagation + relative inference + output propagation",
        "full",
        "propagated",
        "supra_input_relative",
    ),
)

# Every one-component edge in the 2^3 hierarchy-stage factorial, plus supra
# versus single domains. These paired contrasts replace the old ambiguous
# filter-plus-propagation rung.
COMPONENT_EDGES: tuple[tuple[str, str], ...] = (
    ("supra", "single"),
    ("supra_input", "supra"),
    ("supra_relative", "supra"),
    ("supra_output", "supra"),
    ("supra_input_relative", "supra_input"),
    ("supra_input_relative", "supra_relative"),
    ("supra_input_output", "supra_input"),
    ("supra_input_output", "supra_output"),
    ("supra_relative_output", "supra_relative"),
    ("supra_relative_output", "supra_output"),
    ("full", "supra_input_relative"),
    ("full", "supra_input_output"),
    ("full", "supra_relative_output"),
)


# --------------------------------------------------------------------------- #
# Selection-stage accounting (P1: "report counts at every selection stage")     #
# --------------------------------------------------------------------------- #
def selection_stage_counts(
    t0_map: Mapping[str, set],
    t1_exp_map: Mapping[str, set],
    all_architecture_proteins: Iterable[str],
    benchmark_all: Mapping[str, Mapping[str, set]],
    benchmark_with_domains: Mapping[str, Mapping[str, set]],
    ic: Mapping[str, float],
    ic_floors: Sequence[float],
) -> list[dict]:
    """Every filter between "the input files" and "the scored cohort", counted.

    The review's objection is specific and correct: an IC floor does not only
    drop *terms*, it drops *proteins* whose truth becomes empty, so IC >= 0, 2 and
    4 are not necessarily measuring the same proteins and a comparison across
    floors is not paired. This returns one row per selection stage, and for the
    IC stages also ``n_dropped_vs_ic0`` and ``pct_of_ic0`` — how much of the
    unfiltered cohort each floor loses — so that non-comparability is a number in
    the table rather than something a reader has to infer.

    ``benchmark_all`` is the no-knowledge benchmark built **without** the
    has-a-domain restriction; ``benchmark_with_domains`` is the one actually
    scored. Both are ``{aspect: {protein: truth}}``.
    """
    rows: list[dict] = [
        {
            "stage": "t0_annotated_proteins",
            "aspect": "-",
            "min_ic": "-",
            "n_proteins": len(t0_map),
            "n_dropped_vs_ic0": "-",
            "pct_of_ic0": "-",
            "note": "proteins with >=1 GO annotation at t0 under the training "
            "evidence filter (the training universe)",
        },
        {
            "stage": "t1_experimental_proteins",
            "aspect": "-",
            "min_ic": "-",
            "n_proteins": len(t1_exp_map),
            "n_dropped_vs_ic0": "-",
            "pct_of_ic0": "-",
            "note": "proteins with >=1 experimental GO annotation at t1 (the gold standard)",
        },
        {
            "stage": "proteins_with_domains",
            "aspect": "-",
            "min_ic": "-",
            "n_proteins": len(set(all_architecture_proteins)),
            "n_dropped_vs_ic0": "-",
            "pct_of_ic0": "-",
            "note": "proteins with >=1 InterPro domain (predictable at all by a domain method)",
        },
    ]
    for aspect in ("BP", "MF", "CC"):
        rows.append(
            {
                "stage": "no_knowledge_candidates",
                "aspect": aspect,
                "min_ic": "-",
                "n_proteins": len(benchmark_all.get(aspect, {})),
                "n_dropped_vs_ic0": "-",
                "pct_of_ic0": "-",
                "note": "no t0 knowledge in this aspect, gained t1 experimental terms",
            }
        )
        rows.append(
            {
                "stage": "no_knowledge_with_domains",
                "aspect": aspect,
                "min_ic": "-",
                "n_proteins": len(benchmark_with_domains.get(aspect, {})),
                "n_dropped_vs_ic0": "-",
                "pct_of_ic0": "-",
                "note": "the scored cohort before any IC floor",
            }
        )
    for aspect in ("BP", "MF", "CC"):
        base = set(benchmark_with_domains.get(aspect, {}))
        for min_ic in ic_floors:
            kept = set(
                tb.filter_by_ic(benchmark_with_domains.get(aspect, {}), ic, min_ic)
            )
            rows.append(
                {
                    "stage": "cohort_at_ic_floor",
                    "aspect": aspect,
                    "min_ic": f"{min_ic:g}",
                    "n_proteins": len(kept),
                    "n_dropped_vs_ic0": len(base) - len(kept),
                    "pct_of_ic0": (
                        round(100.0 * len(kept) / len(base), 1) if base else 0.0
                    ),
                    "note": "proteins whose truth is non-empty at this IC floor",
                }
            )
    return rows


# --------------------------------------------------------------------------- #
# Score-table helpers                                                          #
# --------------------------------------------------------------------------- #
def rung_prediction_file(
    run_root: Path, rung: Rung, domain_key: str = "interpro"
) -> Path:
    """Where a rung's domain->term table lives.

    The pipeline bakes a non-default domain key into the file stem
    (``domain_ssf_go_...``), mirroring run_dcgo_human's naming.
    """
    prefix = "domain_go" if domain_key == "interpro" else f"domain_{domain_key}_go"
    stem = (
        f"{prefix}_associations_significant.tsv"
        if rung.kind == "associations"
        else f"{prefix}_annotations_propagated.tsv"
    )
    return run_root / rung.run_dir / stem


def run_setting_mismatches(
    run_root: Path, rungs: Iterable[Rung], expected: Mapping[str, str]
) -> list[str]:
    """Rungs whose run manifest disagrees with the evaluator's settings.

    The no-knowledge cohort and IC are built from what training saw, so the
    evaluator's ``domain_key``/``evidence_filter`` must be the ones each run
    used. Rungs without a manifest are not checked.
    """
    problems = []
    for rung in rungs:
        manifest = run_root / rung.run_dir / "run_manifest_go.json"
        if not manifest.exists():
            continue
        params = json.loads(manifest.read_text()).get("parameters", {})
        for key, want in expected.items():
            got = params.get(key)
            if got is not None and got != want:
                problems.append(f"{rung.name}: run {key}={got}, evaluator {want}")
    return problems


def association_count(path: Path) -> int:
    """Rows in a TSV, excluding the header (cheap; the files are large)."""
    with path.open() as handle:
        return max(0, sum(1 for _ in handle) - 1)


# --------------------------------------------------------------------------- #
# Evaluation inputs (shared with check_evaluator_cafaeval.py)                   #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class BenchmarkInputs:
    """Everything the evaluation needs that does not depend on a rung."""

    get_ancestors: Callable[[str], set[str]]
    term_aspect: dict[str, str]
    #: retired GO id -> live id (``OntologyProcessor.alt_id_map`` / ``replaced_by_map``)
    alt_id_map: dict[str, str]
    replaced_by_map: dict[str, str]
    t0_map: dict[str, set[str]]
    t1_exp_map: dict[str, set[str]]
    protein_domains: dict[str, list[str]]
    #: ``{aspect: {protein: truth}}`` without / with the has-a-domain restriction
    benchmark_all: dict[str, dict[str, set[str]]]
    benchmark: dict[str, dict[str, set[str]]]
    eval_proteins: set[str]
    ic: dict[str, float]
    #: propagated t0 term frequency (the naive baseline's scores)
    term_freq: dict[str, float]


def load_benchmark_inputs(
    t0_gaf: Path,
    t1_gaf: Path,
    interpro: Path,
    go_ontology: Path,
    domain_key: str,
    evidence_filter: str,
) -> BenchmarkInputs:  # pragma: no cover - I/O wiring
    """Parse the GAFs, ontology and architectures into the scored benchmark."""
    from loguru import logger

    if str(_ROOT) not in sys.path:
        sys.path.insert(0, str(_ROOT))

    from src.domain_annotation_parser import DomainAnnotationParser
    from src.goa_parser import EXPERIMENTAL_EVIDENCE, GOAParser, parse_goa
    from src.ontology_processor import OntologyProcessor

    logger.info("Loading GO ontology...")
    processor = OntologyProcessor(go_ontology)
    get_ancestors = processor.get_ancestors
    term_aspect = tb.build_term_aspect(processor)

    logger.info(f"Parsing t0 (training) GOA — {evidence_filter} evidence...")
    t0_map = parse_goa(t0_gaf, evidence_filter=evidence_filter)

    logger.info("Parsing t1 (test) GOA — experimental evidence only...")
    t1_exp_map = GOAParser(
        evidence_codes=EXPERIMENTAL_EVIDENCE, aspects={"P", "F", "C"}
    ).parse_gaf_file(t1_gaf)
    t0_map = tb.remap_retired_ids(
        t0_map, processor.alt_id_map, processor.replaced_by_map
    )
    t1_exp_map = tb.remap_retired_ids(
        t1_exp_map, processor.alt_id_map, processor.replaced_by_map
    )

    logger.info("Parsing domain architectures...")
    dom_parser = DomainAnnotationParser(
        max_supra_domain_length=3,
        min_domain_length=10,
        domain_key=domain_key,
    )
    architectures = dom_parser.parse_protein2ipr_file(interpro)
    protein_domains: dict[str, list[str]] = {}
    for protein, arch in architectures.items():
        domains = list(arch.single_domains) + list(arch.supra_domains)
        if domains:
            protein_domains[protein] = domains

    # One cohort for every rung. Supra-domains are built out of single domains,
    # so "has >=1 feature" is the same set either way — which is what makes the
    # paired bootstrap across rungs legitimate.
    logger.info("Building the CAFA no-knowledge benchmark...")
    benchmark_all = tb.build_nk_benchmark_by_aspect(
        t0_map, t1_exp_map, term_aspect, get_ancestors, predictable_proteins=None
    )
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

    logger.info("Computing information content from t0...")
    ic = tb.information_content(t0_map, get_ancestors)
    n_t0 = len(t0_map)
    freq_counts: dict[str, int] = defaultdict(int)
    for terms in t0_map.values():
        for t in tb.propagate_terms(terms, get_ancestors):
            freq_counts[t] += 1
    term_freq = {t: c / n_t0 for t, c in freq_counts.items()} if n_t0 else {}

    return BenchmarkInputs(
        get_ancestors=get_ancestors,
        term_aspect=term_aspect,
        alt_id_map=processor.alt_id_map,
        replaced_by_map=processor.replaced_by_map,
        t0_map=t0_map,
        t1_exp_map=t1_exp_map,
        protein_domains=protein_domains,
        benchmark_all=benchmark_all,
        benchmark=benchmark,
        eval_proteins=eval_proteins,
        ic=ic,
        term_freq=term_freq,
    )


#: ``--transfer`` choice -> domain-to-protein transfer function
TRANSFERS = {
    "pscore": tb.transfer_predictions_pscore,
    "max": tb.transfer_predictions,
}


def build_methods(
    run_root: Path,
    inputs: BenchmarkInputs,
    transfer_name: str = "pscore",
    domain_key: str = "interpro",
) -> tuple[dict[str, dict], list[dict], dict[str, dict]]:  # pragma: no cover - I/O
    """Per-protein predictions for every scored method.

    Each rung under ``-log10(q)``, the ``-log10(p)`` sensitivity variant of each
    association rung (``<rung>__p``), and the naive baseline. Returns
    ``(methods, provenance rows, {rung: domain -> term -log10(q) scores})``;
    the last feeds the permutation null. Raises ``FileNotFoundError`` when a
    rung's predictions are missing.
    """
    transfer = TRANSFERS[transfer_name]
    eval_domains = {p: inputs.protein_domains[p] for p in inputs.eval_proteins}

    def on_scoring_vocabulary(scores: dict) -> dict:
        return tb.remap_retired_scores(
            scores, inputs.alt_id_map, inputs.replaced_by_map, inputs.term_aspect
        )

    methods: dict[str, dict] = {}
    provenance: list[dict] = []
    rung_scores: dict[str, dict] = {}
    for rung in LADDER:
        path = rung_prediction_file(run_root, rung, domain_key)
        if not path.exists():
            raise FileNotFoundError(f"[{rung.name}] missing predictions: {path}")
        n_rows = association_count(path)
        # Primary score: -log10(q). Uniform across rungs, because the True Path
        # output has no p_value column and a ladder scored on two different
        # columns would confound the component with the score scale.
        scores_q = tb.load_domain_go_scores(
            path,
            "q_value" if rung.kind == "propagated" else "adj_p_value",
            neg_log10=True,
        )
        n_domains = len(scores_q)
        scores_q = on_scoring_vocabulary(scores_q)
        rung_scores[rung.name] = scores_q
        methods[rung.name] = transfer(eval_domains, scores_q, inputs.get_ancestors)
        provenance.append(
            {
                "rung": rung.name,
                "label": rung.label,
                "adds": rung.adds,
                "predictions_file": str(path),
                "n_rows": n_rows,
                "n_domains": n_domains,
                "score": "-log10(q)",
            }
        )
        if rung.kind == "associations":
            scores_p = on_scoring_vocabulary(
                tb.load_domain_go_scores(path, "p_value", neg_log10=True)
            )
            methods[f"{rung.name}__p"] = transfer(
                eval_domains, scores_p, inputs.get_ancestors
            )
            provenance.append(
                {
                    "rung": f"{rung.name}__p",
                    "label": f"{rung.label} (p-ranked)",
                    "adds": "sensitivity: score column only",
                    "predictions_file": str(path),
                    "n_rows": n_rows,
                    "n_domains": n_domains,
                    "score": "-log10(p)",
                }
            )
    methods["naive"] = tb.naive_predictions(inputs.eval_proteins, inputs.term_freq)
    provenance.append(
        {
            "rung": "naive",
            "label": "CAFA naive baseline",
            "adds": "-",
            "predictions_file": "-",
            "n_rows": len(inputs.term_freq),
            "n_domains": 0,
            "score": "propagated t0 term frequency",
        }
    )
    return methods, provenance, rung_scores


def floor_cell(
    truth: Mapping[str, set[str]],
    preds: Mapping[str, Mapping[str, float]],
    ic: Mapping[str, float],
    min_ic: float,
) -> tuple[dict, dict, list[float]]:
    """One aspect x IC-floor cell exactly as it is scored.

    The floor restricts **both** sides to terms with ``IC >= min_ic``; proteins
    left with no truth leave the cohort (proteins left with no prediction stay in
    it, as misses). The threshold sweep is drawn from the filtered predictions.
    Returns ``(truth, predictions, thresholds)``.
    """
    pred_f = tb.filter_by_ic(preds, ic, min_ic)
    return tb.filter_by_ic(truth, ic, min_ic), pred_f, tb._candidate_thresholds(pred_f)


def bootstrap_seed(seed: int, aspect: str, min_ic: float) -> int:
    """The paired bootstrap's seed for one aspect x IC-floor cell.

    Derived from ``--seed`` alone, so a re-analysis of the committed cells
    (``check_metric_conventions.py``) can score exactly the resamples their
    intervals were computed on.
    """
    import numpy as np

    rng = np.random.default_rng(seed)
    per_aspect = {a: int(rng.integers(0, 2**31 - 1)) for a in ("BP", "MF", "CC")}
    return per_aspect[aspect] + int(min_ic * 1000)


# --------------------------------------------------------------------------- #
# Driver                                                                       #
# --------------------------------------------------------------------------- #
def main() -> int:  # pragma: no cover - I/O wiring
    import argparse

    import pandas as pd
    from loguru import logger

    logger.remove()
    logger.add(sys.stderr, level="INFO")

    parser = argparse.ArgumentParser(
        description="dcGO component ablation with bootstrap CIs and a permutation null "
        "(VALIDATION_PLAN §4)."
    )
    parser.add_argument("--t0-gaf", type=Path, required=True)
    parser.add_argument("--t1-gaf", type=Path, required=True)
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
        "--domain-key",
        choices=["interpro", "ssf"],
        default="interpro",
        help="Which protein2ipr column defines a domain — must match the "
        "--domain-key the pipeline runs under --run-dir were made with "
        "(default: interpro)",
    )
    parser.add_argument(
        "--evidence-filter",
        choices=["all", "manual", "experimental"],
        default="manual",
        help="t0 evidence filter — must match the --evidence-filter the pipeline "
        "runs were trained with, since the no-knowledge cohort and IC are built "
        "from what training saw (default: manual)",
    )
    parser.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    parser.add_argument("--output-dir", type=Path, default=Path("validation"))
    parser.add_argument("--transfer", choices=sorted(TRANSFERS), default="pscore")
    parser.add_argument(
        "--min-ic",
        type=float,
        action="append",
        metavar="BITS",
        help="Information-content floor(s). Default: 0, 2, 4.",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=1000,
        help="Protein-level bootstrap replicates (default: 1000)",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=200,
        help="Seeded domain-label permutations for the null (default: 200). The "
        "smallest attainable empirical p is 1/(n+1).",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--ci-level", type=float, default=0.95, help="Percentile interval coverage"
    )
    parser.add_argument(
        "--null-rung",
        default="supra",
        help="Which rung the permutation null is computed against (default: supra, "
        "the configuration reported in §2)",
    )
    args = parser.parse_args()
    ic_floors = sorted(set(args.min_ic)) if args.min_ic else [0.0, 2.0, 4.0]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for path in (args.t0_gaf, args.t1_gaf, args.interpro, args.go_ontology):
        if not path.exists():
            logger.error(f"Missing required input: {path}")
            return 1
    mismatches = run_setting_mismatches(
        args.run_dir,
        LADDER,
        {"domain_key": args.domain_key, "evidence_filter": args.evidence_filter},
    )
    if mismatches:
        for problem in mismatches:
            logger.error(f"Settings mismatch — {problem}")
        return 1

    # ---------------------------------------------------------------- inputs --
    inputs = load_benchmark_inputs(
        args.t0_gaf,
        args.t1_gaf,
        args.interpro,
        args.go_ontology,
        args.domain_key,
        args.evidence_filter,
    )
    benchmark, ic, term_aspect = inputs.benchmark, inputs.ic, inputs.term_aspect
    if not inputs.eval_proteins:
        logger.error("Empty benchmark.")
        return 1

    counts = selection_stage_counts(
        inputs.t0_map,
        inputs.t1_exp_map,
        inputs.protein_domains,
        inputs.benchmark_all,
        benchmark,
        ic,
        ic_floors,
    )
    counts_file = args.output_dir / "ablation_selection_counts.tsv"
    pd.DataFrame(counts).to_csv(counts_file, sep="\t", index=False)
    logger.info(f"✓ Selection-stage counts: {counts_file}")

    # ----------------------------------------------------- transfer per rung --
    try:
        methods, provenance, rung_scores = build_methods(
            args.run_dir, inputs, args.transfer, args.domain_key
        )
    except FileNotFoundError as missing:
        logger.error(str(missing))
        return 1
    pd.DataFrame(provenance).to_csv(
        args.output_dir / "ablation_provenance.tsv", sep="\t", index=False
    )

    # --------------------------------------------------- evaluate + bootstrap --
    metric_rows: list[dict] = []
    paired_rows: list[dict] = []

    aspect_preds = {
        name: {
            aspect: tb.restrict_to_aspect(preds, aspect, term_aspect)
            for aspect in ("BP", "MF", "CC")
        }
        for name, preds in methods.items()
    }

    for aspect in ("BP", "MF", "CC"):
        if not benchmark[aspect]:
            continue
        for min_ic in ic_floors:
            true_a = tb.filter_by_ic(benchmark[aspect], ic, min_ic)
            if not true_a:
                continue
            panels: dict[str, rs.EvaluationPanel] = {}
            observed: dict[str, dict[str, float]] = {}
            rows_here: dict[str, dict] = {}
            for name in methods:
                truth_f, pred_f, taus = floor_cell(
                    benchmark[aspect], aspect_preds[name][aspect], ic, min_ic
                )
                panel = rs.build_panel(pred_f, truth_f, ic, taus)
                panels[name] = panel
                observed[name] = rs.panel_metrics(panel)
                row = {
                    "aspect": aspect,
                    "min_ic": min_ic,
                    "method": name,
                    **observed[name],
                }
                rows_here[name] = row
                metric_rows.append(row)
                logger.info(
                    f"  [{aspect} IC>={min_ic:g}] {name:16s} F_max={observed[name]['f_max']:.4f} "
                    f"AUPRC={observed[name]['auprc']:.4f} "
                    f"cov={observed[name]['coverage_at_fmax']:.2f} "
                    f"(n={observed[name]['n_eval_proteins']})"
                )

            logger.info(
                f"  [{aspect} IC>={min_ic:g}] paired bootstrap, {args.n_bootstrap} replicates..."
            )
            reps = rs.paired_bootstrap(
                panels,
                metrics=("f_max", "auprc"),
                n_replicates=args.n_bootstrap,
                seed=bootstrap_seed(args.seed, aspect, min_ic),
                level=args.ci_level,
            )
            for name in panels:
                for metric in ("f_max", "auprc"):
                    lo, hi = rs.percentile_ci(reps[f"{name}::{metric}"], args.ci_level)
                    rows_here[name][f"{metric}_ci_lo"] = lo
                    rows_here[name][f"{metric}_ci_hi"] = hi

            # Paired comparisons that matter: each rung against the component it
            # adds to, and each rung against the naive baseline.
            comparisons = list(COMPONENT_EDGES)
            comparisons += [(r.name, "naive") for r in LADDER]
            comparisons += [("full", "single")]
            for a_name, b_name in comparisons:
                if a_name not in panels or b_name not in panels:
                    continue
                for metric in ("f_max", "auprc"):
                    paired_rows.append(
                        {
                            "aspect": aspect,
                            "min_ic": min_ic,
                            "n_eval_proteins": panels[a_name].n_proteins,
                            **rs.summarise_paired(
                                reps,
                                a_name,
                                b_name,
                                metric,
                                observed[a_name][metric],
                                observed[b_name][metric],
                                level=args.ci_level,
                            ),
                        }
                    )

    metrics_file = args.output_dir / "ablation_metrics.tsv"
    pd.DataFrame(metric_rows).to_csv(metrics_file, sep="\t", index=False)
    logger.info(f"✓ Ablation metrics: {metrics_file}")

    paired_file = args.output_dir / "ablation_paired_bootstrap.tsv"
    pd.DataFrame(paired_rows).to_csv(paired_file, sep="\t", index=False)
    logger.info(f"✓ Paired bootstrap: {paired_file}")

    # ------------------------------------------------------ permutation null --
    if args.n_permutations > 0:
        null_rung = args.null_rung
        if null_rung not in rung_scores:
            logger.error(f"--null-rung {null_rung} is not a rung of the ladder")
            return 1
        logger.info(
            f"Permutation null: {args.n_permutations} seeded domain-label permutations "
            f"of the '{null_rung}' association table..."
        )
        base_scores = rung_scores[null_rung]
        transfer = TRANSFERS[args.transfer]
        eval_domains = {p: inputs.protein_domains[p] for p in inputs.eval_proteins}
        null_samples: dict[tuple[str, float, str], list[float]] = defaultdict(list)
        for i in range(args.n_permutations):
            seed = args.seed * 1_000_003 + i
            shuffled = tb.shuffle_domain_go(base_scores, seed=seed)
            preds = transfer(eval_domains, shuffled, inputs.get_ancestors)
            for aspect in ("BP", "MF", "CC"):
                if not benchmark[aspect]:
                    continue
                pred_a = tb.restrict_to_aspect(preds, aspect, term_aspect)
                for min_ic in ic_floors:
                    true_a = tb.filter_by_ic(benchmark[aspect], ic, min_ic)
                    if not true_a:
                        continue
                    truth_f, pred_f, taus = floor_cell(
                        benchmark[aspect], pred_a, ic, min_ic
                    )
                    panel = rs.build_panel(pred_f, truth_f, ic, taus)
                    vals = rs.panel_metrics(panel)
                    null_samples[(aspect, min_ic, "f_max")].append(vals["f_max"])
                    null_samples[(aspect, min_ic, "auprc")].append(vals["auprc"])
            if (i + 1) % 10 == 0:
                logger.info(f"  permutation {i + 1}/{args.n_permutations}")

        null_rows = []
        observed_lookup = {
            (r["aspect"], r["min_ic"], r["method"]): r for r in metric_rows
        }
        for (aspect, min_ic, metric), samples in sorted(null_samples.items()):
            obs = observed_lookup[(aspect, min_ic, null_rung)][metric]
            null_rows.append(
                {
                    "aspect": aspect,
                    "min_ic": min_ic,
                    "metric": metric,
                    "rung": null_rung,
                    "n_eval_proteins": observed_lookup[(aspect, min_ic, null_rung)][
                        "n_eval_proteins"
                    ],
                    **rs.summarise_null(obs, samples, level=args.ci_level),
                }
            )
        null_file = args.output_dir / "ablation_permutation_null.tsv"
        pd.DataFrame(null_rows).to_csv(null_file, sep="\t", index=False)
        logger.info(f"✓ Permutation null: {null_file}")

    logger.info("=" * 70)
    logger.info("ABLATION COMPLETE")
    for aspect in ("BP", "MF", "CC"):
        for min_ic in ic_floors:
            line = [f"{aspect} IC>={min_ic:g}:"]
            for rung in LADDER:
                row = next(
                    (
                        r
                        for r in metric_rows
                        if r["aspect"] == aspect
                        and r["min_ic"] == min_ic
                        and r["method"] == rung.name
                    ),
                    None,
                )
                if row:
                    line.append(f"{rung.name}={row['f_max']:.3f}")
            if len(line) > 1:
                logger.info("  " + " ".join(line))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
