"""Knowledge-projection benchmark: calibrated unit projection against old dcGO.

Scores the methods of ``src/knowledge_projection.py`` and, when given a §4 run
directory, every rung of the old dcGO ladder, on the CAFA no-knowledge
temporal cohort of ``validation/ablation.py``, with the pre-registered
knowledge-yield criterion (``KNOWLEDGE_PROJECTION.md`` §1.2) and the existing
protein-centric metrics.

Two modes:

- ``temporal`` (default): t0 training, t1 experimental truth, the §4 cohort.
- ``cv``: no t1 at all — t0 aspect-annotated proteins are split into folds and
  each fold is predicted from the others. For choosing the method's own
  settings without touching any held-out label.

    uv run python validation/projection_benchmark.py \\
        --t0-gaf data/raw/goa_archive/goa_human.gaf.205.gz \\
        --t1-gaf data/raw/goa_annotations/goa_human.gaf.gz \\
        --interpro data/interim/protein2ipr_human_t0_ipr85.dat.gz \\
        --cohort-interpro data/interim/protein2ipr_human.dat.gz \\
        --old-run-dir results/ablation-t0ipr-2026-10-15/ipr_manual \\
        --literature data/processed/literature/medembed_protein_y2020.npz \\
        --output-dir validation/projection/final
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
import time
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src import knowledge_projection as kp  # noqa: E402


def _load_sibling(name: str):
    spec = importlib.util.spec_from_file_location(name, _HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault(name, module)
    spec.loader.exec_module(module)
    return module


ab = _load_sibling("ablation")
tb = _load_sibling("temporal_benchmark")
ky = _load_sibling("knowledge_yield")

ASPECTS = ("BP", "MF", "CC")


# --------------------------------------------------------------------------- #
# Unit families                                                                #
# --------------------------------------------------------------------------- #
def unit_families(
    protein_domains: Mapping[str, list[str]], interpro: Path, proteins: set[str]
) -> dict[str, dict[str, set[str]]]:
    """Every unit family, keyed by name, as ``{protein: units}``."""
    singles = {p: {d for d in ds if "," not in d} for p, ds in protein_domains.items()}
    supra = {p: {d for d in ds if "," in d} for p, ds in protein_domains.items()}
    return {
        "ipr": {p: u for p, u in singles.items() if u},
        "supra": {p: u for p, u in supra.items() if u},
        "arch": kp.architecture_units(protein_domains),
        "sig": kp.read_protein_signatures(interpro, proteins),
    }


def combine_families(
    families: Mapping[str, Mapping[str, set[str]]], names: list[str]
) -> dict[str, set[str]]:
    return kp.merge_units(*(families[n] for n in names))


# --------------------------------------------------------------------------- #
# New-method predictions for one aspect                                         #
# --------------------------------------------------------------------------- #
def unit_scores(
    training: kp.AspectTraining,
    units: Mapping[str, set[str]],
    queries: list[str],
    strengths: list[float],
    quantile: float | None = None,
) -> dict[float, sp.csr_matrix]:
    """``{s: query x term max-θ scores}``; the counts are built once."""
    vocab = kp.UnitVocabulary(units)
    evidence = kp.unit_term_evidence(vocab.incidence(training.proteins), training)
    query_incidence = vocab.incidence(queries)
    return {
        s: kp.max_project(query_incidence, evidence.posterior(s, quantile))
        for s in strengths
    }


def literature_scores(
    training: kp.AspectTraining,
    embeddings: Mapping[str, np.ndarray],
    queries: list[str],
    k: int,
    temperature: float,
    max_sim: float = 0.9999,
):
    train_rows = [i for i, p in enumerate(training.proteins) if p in embeddings]
    q_rows = [i for i, p in enumerate(queries) if p in embeddings]
    out = sp.lil_matrix((len(queries), len(training.terms)))
    if not train_rows or not q_rows:
        return out.tocsr()
    train_emb = np.stack([embeddings[training.proteins[i]] for i in train_rows])
    q_emb = np.stack([embeddings[queries[i]] for i in q_rows])
    scores = kp.knn_project(
        q_emb,
        train_emb,
        training.labels[train_rows],
        k=k,
        temperature=temperature,
        max_sim=max_sim,
    )
    expand = sp.csr_matrix(
        (np.ones(len(q_rows)), (q_rows, np.arange(len(q_rows)))),
        shape=(len(queries), len(q_rows)),
    )
    return (expand @ scores).tocsr()


# --------------------------------------------------------------------------- #
# Scoring                                                                      #
# --------------------------------------------------------------------------- #
def leaf_yield(truth, preds, tau: float, get_ancestors) -> float:
    """Correct predictions at ``score ≥ tau`` that have no correct descendant.

    KY counts a correct leaf and every correct ancestor it implies; this
    counts each correct *chain* once, by its most specific member. 0 when
    the precision level was never reached.
    """
    if not np.isfinite(tau):
        return 0.0
    total = 0
    for p, scores in preds.items():
        hits = {t for t, v in scores.items() if v >= tau and t in truth.get(p, ())}
        implied = set()
        for t in hits:
            implied.update(get_ancestors(t))
        total += len(hits - implied)
    return float(total)


def score_methods(
    truth: Mapping[str, set[str]],
    methods: Mapping[str, Mapping[str, Mapping[str, float]]],
    ic: Mapping[str, float],
    min_ic: float,
    protein_metrics: bool,
    get_ancestors=None,
    protein_metrics_for: str | None = None,
    frozen_tau: Mapping[str, float] | None = None,
) -> tuple[list[dict], dict[str, ky.PooledRanking], np.ndarray]:
    """One aspect x IC-floor cell: KY family, pooled AUPRC, CAFA metrics."""
    truth_f = tb.filter_by_ic(truth, ic, min_ic)
    order = sorted(truth_f)
    sizes = np.array([len(truth_f[p]) for p in order], dtype=np.float64)
    rows, rankings = [], {}
    for name, preds in methods.items():
        pred_f = tb.filter_by_ic({p: preds[p] for p in order if p in preds}, ic, min_ic)
        ranking = ky.pool(truth_f, pred_f, ic, proteins=order)
        rankings[name] = ranking
        row = {"method": name, "min_ic": min_ic, "n_proteins": len(order)}
        row |= ky.summarise(ranking)
        if get_ancestors is not None:
            row["ky50_leaf"] = leaf_yield(truth_f, pred_f, row["tau50"], get_ancestors)
        tau = (frozen_tau or {}).get(name)
        if tau is not None and np.isfinite(tau):
            hit = ranking.scores >= tau
            n_pred, n_tp = int(hit.sum()), int(ranking.is_tp[hit].sum())
            row |= {
                "tau_frozen": tau,
                "tp_at_frozen": n_tp,
                "pred_at_frozen": n_pred,
                "precision_at_frozen": n_tp / n_pred if n_pred else float("nan"),
            }
        row["coverage_any"] = sum(1 for p in order if pred_f.get(p)) / max(
            len(order), 1
        )
        if protein_metrics and (
            protein_metrics_for is None or re.search(protein_metrics_for, name)
        ):
            row |= {
                k: v
                for k, v in tb.evaluate_aspect(pred_f, truth_f, ic).items()
                if k in ("f_max", "s_min", "auprc")
            }
        rows.append(row)
    return rows, rankings, sizes


# --------------------------------------------------------------------------- #
# Driver                                                                       #
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--mode", choices=("temporal", "cv"), default="temporal")
    p.add_argument("--t0-gaf", type=Path, required=True)
    p.add_argument("--t1-gaf", type=Path, help="required in temporal mode")
    p.add_argument("--interpro", type=Path, required=True)
    p.add_argument("--cohort-interpro", type=Path)
    p.add_argument(
        "--go-ontology", type=Path, default=Path("data/raw/go_ontology/go-basic.obo")
    )
    p.add_argument("--evidence-filter", default="manual")
    p.add_argument("--old-run-dir", type=Path, help="§4 ladder run dir (old dcGO)")
    p.add_argument("--literature", type=Path, help="literature embeddings .npz")
    p.add_argument(
        "--strength", type=float, nargs="+", default=[1.0], help="Beta prior strength s"
    )
    p.add_argument(
        "--unit-sets",
        nargs="+",
        default=["ipr", "ipr+supra", "ipr+supra+arch+sig"],
        help="unit families joined by '+'",
    )
    p.add_argument(
        "--quantile",
        type=float,
        nargs="+",
        default=[0.0],
        help="posterior lower-bound quantile(s); 0 = posterior mean",
    )
    p.add_argument("--knn-k", type=int, nargs="+", default=[25])
    p.add_argument("--knn-temperature", type=float, nargs="+", default=[0.05])
    p.add_argument(
        "--seq-embeddings", type=Path, nargs="+", help="ProtT5 per-protein .h5 files"
    )
    p.add_argument("--seq-knn-k", type=int, nargs="+", default=[10])
    p.add_argument("--seq-knn-temperature", type=float, nargs="+", default=[0.05])
    p.add_argument(
        "--mix",
        type=lambda v: {
            k: float(x) for k, x in (item.split("=") for item in v.split(","))
        },
        nargs="*",
        default=[],
        help="weighted-mean combinations, e.g. proj=1,lit=0.5 (channels: "
        "proj, lr, seq, lit)",
    )
    p.add_argument(
        "--mix-skip-missing",
        action="store_true",
        help="weighted mean over the channels a protein has data for (sensitivity)",
    )
    p.add_argument(
        "--no-noisy-or",
        action="store_true",
        help="emit only single channels and --mix combinations (the frozen arms)",
    )
    p.add_argument("--lr-units", help="unit families for LR-InterPro, e.g. ipr+sig")
    p.add_argument("--lr-c", type=float, nargs="+", default=[1.0])
    p.add_argument("--lr-min-positives", type=int, default=5)
    p.add_argument("--n-jobs", type=int, default=48)
    p.add_argument("--min-ic", type=float, nargs="+", default=[0.0, 2.0])
    p.add_argument("--aspects", nargs="+", default=list(ASPECTS))
    p.add_argument("--n-folds", type=int, default=5)
    p.add_argument("--n-bootstrap", type=int, default=0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument(
        "--frozen-tau-from",
        type=Path,
        help="projection_metrics.tsv of an earlier split: report precision and "
        "yield at that split's KY50 cutoff per method",
    )
    p.add_argument(
        "--max-proteins-per-paper",
        type=int,
        help="sensitivity: drop t1 truth lines whose papers all annotate more "
        "proteins than this (Schnoes et al. 2013 use 99)",
    )
    p.add_argument(
        "--drop-gene-known",
        action="store_true",
        help="sensitivity: drop cohort proteins whose gene had a same-aspect "
        "manual annotation at t0 under any accession",
    )
    p.add_argument(
        "--save-curves",
        help="regex: write pooled precision/yield curves for matching methods "
        "(plus the best old rung) to projection_curves.tsv",
    )
    p.add_argument(
        "--strict-no-knowledge",
        action="store_true",
        help="sensitivity: keep only cohort proteins with no curated t0 "
        "annotation in any aspect",
    )
    p.add_argument(
        "--reselect-min-ic",
        type=float,
        default=2.0,
        help="IC floors from which the bootstrap re-chooses the best old variant "
        "inside each resample (all 28 variants are resampled)",
    )
    p.add_argument(
        "--cluster-bootstrap",
        action="store_true",
        help="also resample clusters of proteins linked by shared t1 evidence papers",
    )
    p.add_argument("--no-protein-metrics", action="store_true")
    p.add_argument(
        "--protein-metrics-min-ic",
        type=float,
        default=0.0,
        help="compute CAFA protein-centric metrics only at IC floors >= this",
    )
    p.add_argument(
        "--protein-metrics-for",
        help="regex: compute the (slow) CAFA protein-centric metrics only for "
        "matching method names",
    )
    p.add_argument("--stacker", action="store_true", help="add stacked methods")
    p.add_argument("--stack-strength", type=float, default=1.0)
    p.add_argument(
        "--stack-crossfit",
        action="store_true",
        help="2-fold cross-fit on this cohort (development splits only)",
    )
    p.add_argument("--stack-train", type=Path, help="dir of saved dev stack rows")
    p.add_argument("--save-stack-rows", type=Path, help="save this cohort's rows")
    p.add_argument(
        "--extra-train-gaf", type=Path, help="all-species manual GAF (background)"
    )
    p.add_argument("--extra-train-interpro", type=Path, help="its protein2ipr subset")
    p.add_argument(
        "--extra-evidence-filter",
        help="evidence filter for the background only (default: --evidence-filter)",
    )
    p.add_argument(
        "--iea-comparator", action="store_true", help="score t0 IEA annotations"
    )
    p.add_argument(
        "--association-only",
        action="store_true",
        help="score only the association-level endpoint (raw and live-id)",
    )
    p.add_argument("--assoc-quantile", type=float, default=0.1)
    p.add_argument("--assoc-floor", type=float, default=0.3)
    p.add_argument("--interpro2go-t0", type=Path)
    p.add_argument("--interpro2go-t1", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args(argv)


def new_methods_for_aspect(
    args,
    training,
    families,
    embeddings,
    queries,
    get_ancestors,
    seq_embeddings=None,
) -> dict[str, dict[str, dict[str, float]]]:
    """Every configured new-method variant for one aspect's queries.

    Three channels, each a query x term score matrix on one probability-like
    scale: ``proj`` (domain units), ``seq`` (sequence-embedding neighbours)
    and ``lit`` (literature neighbours). Each is reported alone, and every
    combination of one variant per channel is reported under noisy-OR and,
    for each ``--mix`` weight triple, as a weighted mean.
    """
    from loguru import logger

    channels: dict[str, dict[str, sp.csr_matrix]] = {
        "proj": {},
        "lr": {},
        "seq": {},
        "lit": {},
    }
    for unit_set in args.unit_sets:
        units = combine_families(families, unit_set.split("+"))
        t = time.time()
        for q in args.quantile:
            qq = None if q <= 0 else q
            qtag = "" if qq is None else f",q={q:g}"
            for s, mat in unit_scores(
                training, units, queries, args.strength, qq
            ).items():
                channels["proj"][f"proj[{unit_set}|s={s:g}{qtag}]"] = mat
        logger.info(f"    proj[{unit_set}]: {time.time() - t:.1f}s")
    if args.lr_units:
        t = time.time()
        vocab = kp.UnitVocabulary(combine_families(families, args.lr_units.split("+")))
        ancestors = kp.ancestor_matrix(training.terms, get_ancestors)
        for C in args.lr_c:
            channels["lr"][f"lr[{args.lr_units}|C={C:g}]"] = kp.lr_project(
                vocab.incidence(training.proteins),
                training,
                vocab.incidence(queries),
                min_positives=args.lr_min_positives,
                C=C,
                n_jobs=args.n_jobs,
                ancestors=ancestors,
            )
        logger.info(f"    lr[{args.lr_units}]: {time.time() - t:.1f}s")
    for name, emb in (("lit", embeddings), ("seq", seq_embeddings)):
        if not emb:
            continue
        ks = args.knn_k if name == "lit" else args.seq_knn_k
        temps = args.knn_temperature if name == "lit" else args.seq_knn_temperature
        for k in ks:
            for temp in temps:
                t = time.time()
                # Identical literature vectors are one gene under two
                # accessions; identical sequences are real (e.g. histone
                # genes), so only literature drops them.
                max_sim = 0.9999 if name == "lit" else 1.01
                mat = literature_scores(training, emb, queries, k, temp, max_sim)
                channels[name][f"{name}[k={k},T={temp:g}]"] = mat
                logger.info(f"    {name}[k={k},T={temp:g}] in {time.time() - t:.1f}s")

    out = {
        n: kp.to_prediction_dict(m, queries, training.terms)
        for ch in channels.values()
        for n, m in ch.items()
    }
    present = [c for c in ("proj", "lr", "seq", "lit") if channels[c]]
    from itertools import combinations, product

    for size in range(2, len(present) + 1):
        for group in combinations(present, size):
            for picks in product(*(channels[c].items() for c in group)):
                names = [n for n, _ in picks]
                mats = [m for _, m in picks]
                both = mats[0]
                for m in mats[1:]:
                    both = kp.noisy_or(both, m)
                if not args.no_noisy_or:
                    out["+".join(names)] = kp.to_prediction_dict(
                        both, queries, training.terms
                    )
                for spec in args.mix:
                    if set(spec) != set(group):
                        continue
                    ws = [spec[c] for c in group]
                    tag = ",".join(f"{c}={spec[c]:g}" for c in group)
                    out[f"mix({tag})[" + "+".join(names) + "]"] = kp.to_prediction_dict(
                        kp.weighted_mean(mats, ws, args.mix_skip_missing),
                        queries,
                        training.terms,
                    )
    return out


def old_methods(run_dir: Path, inputs) -> dict[str, dict]:
    """Every old-ladder rung under both transfers and both score columns."""
    out: dict[str, dict] = {}
    eval_domains = {p: inputs.protein_domains.get(p, []) for p in inputs.eval_proteins}
    for rung in ab.LADDER:
        path = ab.rung_prediction_file(run_dir, rung, "interpro")
        columns = ["q_value" if rung.kind == "propagated" else "adj_p_value"]
        if rung.kind == "associations":
            columns.append("p_value")
        for col in columns:
            scores = tb.remap_retired_scores(
                tb.load_domain_go_scores(path, col, neg_log10=True),
                inputs.alt_id_map,
                inputs.replaced_by_map,
                inputs.term_aspect,
            )
            tag = "q" if col != "p_value" else "p"
            for tname, transfer in ab.TRANSFERS.items():
                out[f"old[{rung.name}|{tag}|{tname}]"] = transfer(
                    eval_domains, scores, inputs.get_ancestors
                )
    return out


def without_high_throughput(args, inputs):
    """The same cohort scored on t1 truth minus high-throughput papers.

    A sensitivity analysis (Schnoes et al. 2013): experimental t1 lines whose
    every cited PMID annotates more than ``--max-proteins-per-paper``
    proteins are dropped before the no-knowledge benchmark is rebuilt.
    Proteins left without truth leave the cohort; nothing else changes.
    """
    import dataclasses
    import gzip

    from loguru import logger

    from src.goa_parser import EXPERIMENTAL_EVIDENCE, has_not_qualifier

    lines = []
    proteins_per_paper: dict[str, set[str]] = {}
    with gzip.open(args.t1_gaf, "rt") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if (
                len(f) < 9
                or f[6] not in EXPERIMENTAL_EVIDENCE
                or has_not_qualifier(f[3])
            ):
                continue
            pmids = [r for r in f[5].split("|") if r.startswith("PMID:")]
            lines.append((f[1], f[4], pmids))
            for r in pmids:
                proteins_per_paper.setdefault(r, set()).add(f[1])
    big = {
        r
        for r, ps in proteins_per_paper.items()
        if len(ps) > args.max_proteins_per_paper
    }
    t1: dict[str, set[str]] = {}
    for protein, term, pmids in lines:
        if pmids and all(r in big for r in pmids):
            continue
        t1.setdefault(protein, set()).add(term)
    t1 = tb.remap_retired_ids(t1, inputs.alt_id_map, inputs.replaced_by_map)
    benchmark = tb.build_nk_benchmark_by_aspect(
        inputs.t0_map,
        t1,
        inputs.term_aspect,
        inputs.get_ancestors,
        predictable_proteins=inputs.cohort_proteins,
    )
    logger.info(
        f"High-throughput filter: {len(big):,} papers annotate > "
        f"{args.max_proteins_per_paper} proteins; cohort "
        + ", ".join(f"{a} {len(benchmark[a]):,}" for a in ASPECTS)
    )
    return dataclasses.replace(
        inputs,
        t1_exp_map=t1,
        benchmark=benchmark,
        eval_proteins={p for a in benchmark.values() for p in a},
    )


def load_frozen_taus(path: Path | None) -> dict[tuple[str, float], dict[str, float]]:
    """``{(aspect, min_ic): {method: tau50}}`` from an earlier split's metrics.

    The deployable reading of KY50: a curator fixes the score cutoff where
    the development split reached 50% precision, then sees what precision
    and yield that cutoff gives on a later split.
    """
    if path is None:
        return {}
    frame = pd.read_csv(path, sep="\t")
    out: dict[tuple[str, float], dict[str, float]] = {}
    for (aspect, min_ic), sub in frame.groupby(["aspect", "min_ic"]):
        out[(aspect, float(min_ic))] = dict(zip(sub["method"], sub["tau50"]))
    return out


def comparators(args, inputs) -> dict[str, dict]:
    """Knowledge already on record at t0, as binary predictions (score 1).

    ``interpro2go_t0``: the t0 InterPro2GO mappings of each protein's entries,
    propagated. ``iea_t0``: every electronic (IEA) annotation the t0 GAF
    already held for the protein, propagated — InterPro2GO, UniRule, ARBA,
    Ensembl Compara and the rest. Neither is used by the new method; they
    say how much of the future was already predicted at t0.
    """
    from src.goa_parser import GOAParser

    out: dict[str, dict] = {}
    roots = kp.ASPECT_ROOTS

    def closure(terms):
        closed = set()
        for t in terms:
            closed.add(t)
            closed.update(inputs.get_ancestors(t))
        return {t: 1.0 for t in closed - roots}

    if args.interpro2go_t0:
        by_entry: dict[str, set[str]] = {}
        for entry, term in ab.tdc.read_interpro2go(args.interpro2go_t0):
            by_entry.setdefault(entry, set()).add(term)
        preds = {}
        for p in inputs.eval_proteins:
            terms = set()
            for d in inputs.protein_domains.get(p, []):
                terms |= by_entry.get(d, set())
            if terms:
                preds[p] = closure(terms)
        out["comparator[interpro2go_t0]"] = preds
    if args.iea_comparator:
        iea = GOAParser(evidence_codes={"IEA"}, aspects={"P", "F", "C"}).parse_gaf_file(
            args.t0_gaf
        )
        iea = tb.remap_retired_ids(iea, inputs.alt_id_map, inputs.replaced_by_map)
        out["comparator[iea_t0]"] = {
            p: closure(iea[p]) for p in inputs.eval_proteins if iea.get(p)
        }
    return out


GAF_ASPECT = {"P": "BP", "F": "MF", "C": "CC"}


def gaf_rows(path: Path):
    """``(accession, symbol, taxon, evidence, aspect, negated)`` per GAF line."""
    import gzip

    from src.goa_parser import has_not_qualifier

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            yield (
                f[1],
                f[2],
                f[12].split("|")[0],
                f[6],
                GAF_ASPECT.get(f[8]),
                has_not_qualifier(f[3]),
            )


def cohort_symbols(args, inputs) -> dict[str, str]:
    """Gene symbol of every cohort protein, from the human t1 and t0 GAFs."""
    out: dict[str, str] = {}
    for path in (args.t1_gaf, args.t0_gaf):
        for acc, symbol, *_ in gaf_rows(path):
            if acc in inputs.eval_proteins:
                out.setdefault(acc, symbol)
    return out


def strict_no_knowledge(inputs):
    """Sensitivity cohort: proteins with no curated t0 annotation in *any* aspect.

    The primary cohort is no-knowledge in the target aspect only (CAFA would
    call proteins annotated in another aspect "limited knowledge"). This
    keeps only proteins the t0 training map does not mention at all.
    """
    import dataclasses

    from loguru import logger

    known = {p for p, terms in inputs.t0_map.items() if terms}
    benchmark = {
        a: {p: t for p, t in inputs.benchmark[a].items() if p not in known}
        for a in ASPECTS
    }
    logger.info(
        "Strict no-knowledge cohort: "
        + ", ".join(f"{a} {len(benchmark[a]):,}" for a in ASPECTS)
    )
    return dataclasses.replace(
        inputs,
        benchmark=benchmark,
        eval_proteins={p for a in benchmark.values() for p in a},
    )


def without_gene_known(args, inputs):
    """Sensitivity cohort: drop proteins whose *gene* was already known at t0.

    A cohort protein has no manual annotation in the aspect under its own
    accession at t0, but its gene may have had one under another accession
    (reference-proteome switches such as HLA-A P01891 → P04439) or in the
    all-species file the human GAF's redundancy filter thins. Predicting
    those is rediscovery, not projection. Per aspect, proteins whose symbol
    has a same-aspect manual (non-IEA, non-NOT) human row in the t0 human
    GAF or ``--extra-train-gaf`` are dropped.
    """
    import dataclasses

    from loguru import logger

    from src.goa_parser import NON_IEA_EVIDENCE

    symbols = cohort_symbols(args, inputs)
    known: dict[str, set[str]] = {a: set() for a in ASPECTS}
    sources = [args.t0_gaf] + ([args.extra_train_gaf] if args.extra_train_gaf else [])
    for path in sources:
        for _acc, symbol, taxon, evidence, aspect, negated in gaf_rows(path):
            if (
                aspect
                and taxon == "taxon:9606"
                and evidence in NON_IEA_EVIDENCE
                and not negated
            ):
                known[aspect].add(symbol)
    benchmark = {
        a: {
            p: t
            for p, t in inputs.benchmark[a].items()
            if symbols.get(p) not in known[a]
        }
        for a in ASPECTS
    }
    logger.info(
        "Gene-known filter: dropped "
        + ", ".join(
            f"{a} {len(inputs.benchmark[a]) - len(benchmark[a]):,}" for a in ASPECTS
        )
    )
    return dataclasses.replace(
        inputs,
        benchmark=benchmark,
        eval_proteins={p for a in benchmark.values() for p in a},
    )


def gaf_accessions(path: Path) -> set[str]:
    """Every accession a GAF lists, under any evidence."""
    import gzip

    out: set[str] = set()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!"):
                out.add(line.split("\t", 2)[1])
    return out


def check_old_runs(args) -> list[str]:
    """The old ladder must have been trained on this split's t0 inputs.

    Same guard as ``ablation.py``: every rung manifest's domain key,
    evidence filter, protein2ipr and t0 GAF digests must match.
    """
    from src.run_manifest import sha256_file

    return ab.run_setting_mismatches(
        args.old_run_dir,
        ab.LADDER,
        {"domain_key": "interpro", "evidence_filter": args.evidence_filter},
        input_sha256={
            "domain_annotations": sha256_file(args.interpro),
            "gaf": sha256_file(args.t0_gaf),
        },
    )


def run_temporal(args) -> int:
    from loguru import logger

    if args.old_run_dir:
        problems = check_old_runs(args)
        for problem in problems:
            logger.error(f"Old ladder mismatch — {problem}")
        if problems:
            return 1

    inputs = ab.load_benchmark_inputs(
        args.t0_gaf,
        args.t1_gaf,
        args.interpro,
        args.go_ontology,
        "interpro",
        args.evidence_filter,
        cohort_interpro=args.cohort_interpro,
    )
    if args.max_proteins_per_paper:
        inputs = without_high_throughput(args, inputs)
    if args.drop_gene_known:
        inputs = without_gene_known(args, inputs)
    if args.strict_no_knowledge:
        inputs = strict_no_knowledge(inputs)
    families = unit_families(
        inputs.protein_domains, args.interpro, set(inputs.protein_domains)
    )
    embeddings, n_papers = (
        kp.load_literature(args.literature) if args.literature else ({}, {})
    )
    logger.info(
        "Unit families: "
        + ", ".join(f"{k}={len(v):,} proteins" for k, v in families.items())
        + (f"; literature {len(embeddings):,} proteins" if embeddings else "")
    )
    old = old_methods(args.old_run_dir, inputs) if args.old_run_dir else {}
    old |= comparators(args, inputs)
    naive = tb.naive_predictions(inputs.eval_proteins, inputs.term_freq)
    train_map, tag = inputs.t0_map, ""
    if args.extra_train_gaf:
        train_map = add_background(args, inputs, families)
        tag = f"@allsp-{args.extra_evidence_filter or args.evidence_filter}"
    seq_embeddings = (
        kp.load_h5_embeddings(
            args.seq_embeddings,
            accessions=set(train_map) | inputs.eval_proteins,
        )
        if args.seq_embeddings
        else {}
    )
    if seq_embeddings:
        logger.info(f"Sequence embeddings: {len(seq_embeddings):,} proteins")
    # Embedding files are current; a cohort protein with no record at t0 (no
    # t0 InterPro match and no line in the t0 GAF) may not have existed then,
    # and old dcGO cannot score it either, so no channel may.
    existed = set(inputs.protein_domains) | gaf_accessions(args.t0_gaf)
    late = inputs.eval_proteins - existed
    embeddings = {p: v for p, v in embeddings.items() if p not in late}
    seq_embeddings = {p: v for p, v in seq_embeddings.items() if p not in late}
    logger.info(f"Cohort proteins with no t0 record (no embedding used): {len(late):,}")

    frozen_taus = load_frozen_taus(args.frozen_tau_from)
    # Queries never train, whatever the background file holds about them.
    trainings = kp.build_aspect_trainings(
        train_map,
        inputs.get_ancestors,
        inputs.term_aspect,
        args.aspects,
        exclude={a: list(inputs.benchmark[a]) for a in args.aspects},
    )
    all_rows, boot_rows, curve_rows = [], [], []
    for aspect in [] if args.association_only else args.aspects:
        truth = inputs.benchmark[aspect]
        if not truth:
            continue
        logger.info(f"[{aspect}] {len(truth):,} cohort proteins")
        queries = sorted(truth)
        training = trainings[aspect]
        methods = new_methods_for_aspect(
            args,
            training,
            families,
            embeddings,
            queries,
            inputs.get_ancestors,
            seq_embeddings,
        )
        if args.stacker:
            methods |= stacker_methods(
                args,
                aspect,
                training,
                families,
                embeddings,
                n_papers,
                queries,
                truth,
                inputs.ic,
                inputs.get_ancestors,
            )
        if tag:
            methods = {f"{name}{tag}": preds for name, preds in methods.items()}
        for name, preds in {**old, "naive": naive}.items():
            methods[name] = tb.restrict_to_aspect(
                {p: preds[p] for p in queries if p in preds}, aspect, inputs.term_aspect
            )
        for min_ic in args.min_ic:
            rows, rankings, sizes = score_methods(
                truth,
                methods,
                inputs.ic,
                min_ic,
                not args.no_protein_metrics and min_ic >= args.protein_metrics_min_ic,
                inputs.get_ancestors,
                args.protein_metrics_for,
                frozen_taus.get((aspect, min_ic)),
            )
            for row in rows:
                row["aspect"] = aspect
            all_rows.extend(rows)
            if args.save_curves:
                curve_rows.extend(
                    curve_points(rankings, rows, aspect, min_ic, args.save_curves)
                )
            if args.n_bootstrap:
                order = next(iter(rankings.values())).proteins
                clusters = (
                    evidence_clusters(args.t1_gaf, aspect, order)
                    if args.cluster_bootstrap and min_ic >= args.reselect_min_ic
                    else None
                )
                boot_rows.extend(
                    bootstrap_rows(
                        rows,
                        rankings,
                        sizes,
                        aspect,
                        min_ic,
                        args,
                        frozen_taus.get((aspect, min_ic)),
                        clusters,
                    )
                )
            logger.info(f"  IC>={min_ic:g}: scored {len(rows)} methods")
    if not args.association_only:
        write_outputs(args, all_rows, boot_rows)
    if curve_rows:
        pd.DataFrame(curve_rows).to_csv(
            args.output_dir / "projection_curves.tsv", sep="\t", index=False
        )
    if args.interpro2go_t0 and args.old_run_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for remap, suffix in ((False, ""), (True, "_live_ids")):
            assoc = pd.DataFrame(association_level(args, inputs, families, remap))
            assoc.to_csv(
                args.output_dir / f"projection_association_level{suffix}.tsv",
                sep="\t",
                index=False,
            )
    return 0


def stacker_methods(
    args,
    aspect,
    training,
    families,
    embeddings,
    n_papers,
    queries,
    truth,
    ic,
    get_ancestors,
) -> dict[str, dict]:
    """Stacked predictions: cross-fitted on this cohort and/or trained on dev rows.

    Two arms, each with its own candidates: ``domain`` sees no literature at
    all (not even through which terms become candidates); ``all`` adds it.
    Rows are saved per arm with ``--save-stack-rows`` so a later split can
    train on them (``--stack-train``).
    """
    from loguru import logger

    from src import projection_stacker as ps

    ancestors = kp.ancestor_matrix(training.terms, get_ancestors)
    arms = {"domain": None}
    if embeddings:
        arms["all"] = embeddings
    use = {
        "domain": [
            f for f in ps.FEATURES if not f.startswith("lit_") and f != "log_n_papers"
        ],
        "all": list(ps.FEATURES),
    }
    out: dict[str, dict] = {}
    for arm, emb in arms.items():
        t = time.time()
        feats = ps.build_features(
            training,
            families,
            queries,
            args.stack_strength,
            ic,
            emb,
            n_papers,
            knn_k=args.knn_k[0],
        )
        y = feats.labels(truth)
        logger.info(
            f"    stack[{arm}]: {feats.rows.size:,} candidates, {int(y.sum()):,} "
            f"positive, built in {time.time() - t:.1f}s"
        )
        if args.save_stack_rows:
            args.save_stack_rows.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(
                args.save_stack_rows / f"stack_rows_{arm}_{aspect}.npz",
                X=feats.X,
                y=y,
                features=np.array(ps.FEATURES),
                rows=feats.rows,
                cols=feats.cols,
                queries=np.array(feats.queries),
                terms=np.array(feats.terms),
            )
        if args.stack_crossfit:
            p = ps.cross_fit(feats, y, use[arm], n_folds=2, seed=args.seed)
            out[f"stack[{arm}|2fold]"] = kp.to_prediction_dict(
                ps.finish(feats, p, ancestors), queries, training.terms
            )
        if args.stack_train:
            rows = np.load(args.stack_train / f"stack_rows_{arm}_{aspect}.npz")
            if tuple(rows["features"]) != ps.FEATURES:
                raise SystemExit("--stack-train rows were built with other features")
            model = ps.fit(rows["X"], rows["y"], use[arm], seed=args.seed)
            out[f"stack[{arm}|dev-trained]"] = kp.to_prediction_dict(
                ps.predict(model, feats, ancestors), queries, training.terms
            )
        logger.info(f"    stack[{arm}] done in {time.time() - t:.1f}s")
    return out


def add_background(args, inputs, families) -> dict[str, set[str]]:
    """Merge an all-species manual GAF into training; add its proteins' units.

    Only units some cohort protein carries are read for background proteins
    (``read_restricted_units``) — nothing else can be projected. Human
    proteins keep the units of the human subset; their terms are the union of
    both files'.
    """
    from loguru import logger

    from src.goa_parser import parse_goa

    extra = tb.remap_retired_ids(
        parse_goa(
            args.extra_train_gaf,
            evidence_filter=args.extra_evidence_filter or args.evidence_filter,
        ),
        inputs.alt_id_map,
        inputs.replaced_by_map,
    )
    cohort = {p for truths in inputs.benchmark.values() for p in truths}
    vocabulary = {
        u for fam in families.values() for p in cohort for u in fam.get(p, ())
    }
    # A cohort protein's own gene must never train for it: drop human
    # background accessions that share a cohort protein's gene symbol.
    cohort_genes = set(cohort_symbols(args, inputs).values())
    same_gene = {
        acc
        for acc, symbol, taxon, *_ in gaf_rows(args.extra_train_gaf)
        if taxon == "taxon:9606" and symbol in cohort_genes and acc not in cohort
    }
    extra = {p: t for p, t in extra.items() if p not in same_gene}
    background = set(extra) - set(inputs.protein_domains)
    t = time.time()
    units = kp.read_restricted_units(args.extra_train_interpro, background, vocabulary)
    for fam, by_protein in kp.split_families(units).items():
        families[fam] = {**by_protein, **families[fam]}
    merged = {p: set(ts) for p, ts in extra.items()}
    for p, ts in inputs.t0_map.items():
        merged.setdefault(p, set()).update(ts)
    logger.info(
        f"Background: {len(extra):,} annotated proteins ({len(same_gene):,} "
        f"same-gene human accessions dropped), {len(units):,} with a cohort unit "
        f"({len(vocabulary):,} units) in {time.time() - t:.0f}s"
    )
    return merged


def confident_leaves(
    evidence, strength, training, vocab, get_ancestors, quantile, floor
) -> list[tuple[str, str, float]]:
    """Per unit, its most specific terms whose lower credible bound clears ``floor``.

    A curated InterPro2GO mapping is a claim that (nearly) every member has
    the term, so the bound is the right score; the evaluator re-adds every
    ancestor itself, so only each confident chain's most specific term is
    listed. Ranked (ascending) by ``1 − bound``.
    """
    lcb = evidence.posterior(strength, quantile).tocsr()
    lcb.data[lcb.data < floor] = 0.0
    lcb.eliminate_zeros()
    strict = kp.ancestor_matrix(training.terms, get_ancestors)
    strict.setdiag(0)
    strict.eliminate_zeros()
    present = (lcb > 0).astype(np.float64)
    implied = (present @ strict) > 0
    leaves = lcb.multiply(present - present.multiply(implied)).tocoo()
    return [
        (vocab.units[u], training.terms[t], 1.0 - v)
        for u, t, v in zip(leaves.row, leaves.col, leaves.data)
        if v > 0
    ]


def association_level(args, inputs, families, remap: bool = False) -> list[dict]:
    """The project's pre-registered primary endpoint, with the new estimator.

    Single InterPro entries only (as pre-registered). The new method's table
    is every (entry, term) with ``k ≥ 2`` annotated carriers, ranked by
    ``θ`` (fed to the evaluator as ``1 − θ``, which ranks ascending). The
    budget is the smallest **old** rung's eligible count, so adding a table
    cannot move it.
    """
    from loguru import logger

    tdc = ab.tdc
    t0_pairs = tdc.read_interpro2go(args.interpro2go_t0)
    t1_pairs = tdc.read_interpro2go(args.interpro2go_t1)
    if remap:
        # A t0 mapping on a GO id retired since, and the same mapping on its
        # successor at t1, is a rename, not an addition: put both snapshots on
        # live ids first (alt_id, then replaced_by — the input cleanup rule).
        def live(pairs):
            out = set()
            for d, t in pairs:
                t = inputs.alt_id_map.get(t, t)
                out.add((d, inputs.replaced_by_map.get(t, t)))
            return out

        t0_pairs, t1_pairs = live(t0_pairs), live(t1_pairs)
    delta = tdc.curated_delta(t0_pairs, t1_pairs, inputs.get_ancestors)
    universe, _ = tdc.read_protein2ipr_universe(args.interpro, "interpro")
    old = [
        tdc.load_rung(r.name, ab.rung_prediction_file(args.old_run_dir, r, "interpro"))
        for r in ab.LADDER
    ]
    budget = tdc.default_budget(old, tdc.eligible_domains(delta, universe))
    vocab = kp.UnitVocabulary(families["ipr"])
    trainings = {
        aspect: kp.build_aspect_training(
            inputs.t0_map, inputs.get_ancestors, inputs.term_aspect, aspect
        )
        for aspect in ASPECTS
    }
    for s in args.strength:
        mean_rows, leaf_rows = [], []
        for aspect, training in trainings.items():
            evidence = kp.unit_term_evidence(
                vocab.incidence(training.proteins), training
            )
            theta = evidence.posterior(s).tocoo()
            k = evidence.counts.tocoo()
            keep = k.data >= 2  # same pattern and order as theta (both from counts)
            mean_rows.extend(
                (vocab.units[u], training.terms[t], 1.0 - v)
                for u, t, v in zip(theta.row[keep], theta.col[keep], theta.data[keep])
            )
            leaf_rows.extend(
                confident_leaves(
                    evidence,
                    s,
                    training,
                    vocab,
                    inputs.get_ancestors,
                    args.assoc_quantile,
                    args.assoc_floor,
                )
            )
        old.append(tdc.RungTable(f"proj[ipr|s={s:g}]", tuple(mean_rows), 0))
        old.append(
            tdc.RungTable(
                f"proj[ipr|s={s:g},q={args.assoc_quantile:g},leaf]",
                tuple(leaf_rows),
                0,
            )
        )
    logger.info(f"Association level: budget {budget:,} (smallest old rung)")
    out = tdc.evaluate_rungs(
        old,
        delta,
        inputs.get_ancestors,
        universe=universe,
        budget=budget,
        reference="supra_input_relative",
        n_bootstrap=max(args.n_bootstrap, 200),
        seed=args.seed,
    )
    return out


def run_cv(args) -> int:
    """Fold-wise prediction of t0 labels; no held-out time point is read."""
    from loguru import logger

    from src.goa_parser import parse_goa
    from src.ontology_processor import OntologyProcessor

    processor = OntologyProcessor(args.go_ontology)
    term_aspect = tb.build_term_aspect(processor)
    t0_map = tb.remap_retired_ids(
        parse_goa(args.t0_gaf, evidence_filter=args.evidence_filter),
        processor.alt_id_map,
        processor.replaced_by_map,
    )
    ic = tb.information_content(t0_map, processor.get_ancestors)
    protein_domains = ab_feature_map(args.interpro)
    families = unit_families(protein_domains, args.interpro, set(protein_domains))
    embeddings = (
        kp.load_literature_embeddings(args.literature) if args.literature else {}
    )
    rng = np.random.default_rng(args.seed)

    all_rows = []
    for aspect in args.aspects:
        full = kp.build_aspect_training(
            t0_map, processor.get_ancestors, term_aspect, aspect
        )
        with_units = [p for p in full.proteins if p in protein_domains]
        folds = rng.integers(0, args.n_folds, size=len(with_units))
        fold_methods: dict[str, dict] = {}
        truth: dict[str, set[str]] = {}
        for f in range(args.n_folds):
            held = [p for p, k in zip(with_units, folds) if k == f]
            training = kp.build_aspect_training(
                t0_map, processor.get_ancestors, term_aspect, aspect, exclude=held
            )
            logger.info(f"[{aspect}] fold {f}: {len(held):,} held out")
            preds = new_methods_for_aspect(
                args, training, families, embeddings, held, processor.get_ancestors
            )
            for name, by_protein in preds.items():
                fold_methods.setdefault(name, {}).update(by_protein)
        labels = full.labels.tocsr()
        index = {p: i for i, p in enumerate(full.proteins)}
        for p in with_units:
            r = labels[index[p]]
            truth[p] = {full.terms[j] for j in r.indices}
        for min_ic in args.min_ic:
            rows, _, _ = score_methods(
                truth, fold_methods, ic, min_ic, not args.no_protein_metrics
            )
            for row in rows:
                row["aspect"] = aspect
            all_rows.extend(rows)
    write_outputs(args, all_rows, [])
    return 0


def ab_feature_map(interpro: Path) -> dict[str, list[str]]:
    from src.domain_annotation_parser import DomainAnnotationParser

    architectures = DomainAnnotationParser(
        max_supra_domain_length=3, min_domain_length=10, domain_key="interpro"
    ).parse_protein2ipr_file(interpro)
    return {
        p: list(a.single_domains) + list(a.supra_domains)
        for p, a in architectures.items()
        if a.single_domains
    }


def curve_points(rankings, rows, aspect, min_ic, pattern, max_points=400):
    """Pooled (yield, precision) at distinct thresholds, for plotting.

    Methods whose name matches ``pattern``, plus the best old rung (the KY50
    reference, chosen as in :func:`bootstrap_rows`). Long curves are thinned
    to ``max_points`` evenly spaced thresholds, always keeping the last.
    """
    old_rows = [r for r in rows if r["method"].startswith("old[")]
    keep = {r["method"] for r in rows if re.search(pattern, r["method"])}
    if old_rows:
        keep.add(
            max(old_rows, key=lambda r: (r["ky50"], r["pooled_auprc"], r["method"]))[
                "method"
            ]
        )
    out = []
    for name in sorted(keep):
        r = rankings[name]
        if r.scores.size == 0:
            continue
        tp, pp, _ = ky._curves(r, None)
        # Every threshold while precision is still >= 20% (where KY30, KY50
        # and KY70 are decided); the long low-precision tail is thinned.
        prec = np.where(pp > 0, tp / np.maximum(pp, 1e-300), 0.0)
        above = np.flatnonzero(prec >= 0.2)
        head = np.arange(above.max() + 1) if above.size else np.arange(0)
        tail = np.linspace(0, tp.size - 1, min(max_points, tp.size)).astype(int)
        idx = np.unique(np.r_[head, tail])
        for i in idx:
            out.append(
                {
                    "aspect": aspect,
                    "min_ic": min_ic,
                    "method": name,
                    "threshold": float(r.scores[r.group_end[i] - 1]),
                    "tp": float(tp[i]),
                    "predicted": float(pp[i]),
                    "precision": float(tp[i] / pp[i]) if pp[i] else 0.0,
                    "n_truth": r.n_truth,
                }
            )
    return out


BOOT_METRICS = ("ky30", "ky50", "ky70", "pooled_auprc")


def best_old_by_metric(rows) -> dict[str, str]:
    """The old rung that does best on each metric, among every old variant.

    KY50's is the pre-registered reference. For the others, comparing against
    the rung chosen by KY50 would flatter the new method, so each metric gets
    the old variant that is best on that metric. Ties are broken by pooled
    AUPRC (KY50 for AUPRC itself), then by name, so the choice is stable.
    """
    old_rows = [r for r in rows if r["method"].startswith("old[")]
    if not old_rows:
        return {}
    out = {}
    for metric in BOOT_METRICS:
        second = "ky50" if metric == "pooled_auprc" else "pooled_auprc"
        out[metric] = max(old_rows, key=lambda r: (r[metric], r[second], r["method"]))[
            "method"
        ]
    return out


def _ratio(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    """Per-replicate ratio: x/0 = inf, 0/0 dropped (as ``ky.ratio_interval``)."""
    keep = (num > 0) | (den > 0)
    num, den = num[keep], den[keep]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, num / np.where(den > 0, den, 1.0), np.inf)


def bootstrap_rows(
    rows,
    rankings,
    sizes,
    aspect,
    min_ic,
    args,
    frozen_tau=None,
    clusters=None,
) -> list[dict]:
    """Paired bootstrap of every new method against old dcGO, metric by metric.

    Each metric's ``reference`` is the old variant best on that metric on the
    full cohort (:func:`best_old_by_metric`). At IC floors from
    ``--reselect-min-ic`` the best old variant is also re-chosen inside every
    resample (``ratio_reselect_*``, ``p_ratio_reselect_ge10``), which carries
    the selection's own uncertainty (``ratio_reselect_*``, ``diff_reselect_*``).
    Differences new − reference, ratios to
    naive (AUPRC), and the precision at a frozen cutoff are reported too.
    With ``clusters`` the whole analysis is repeated resampling clusters of
    proteins instead of proteins (``resampling`` column).
    """
    refs = best_old_by_metric(rows)
    if not refs:
        return []
    olds = [n for n in rankings if n.startswith("old[")]
    news = [n for n in rankings if not n.startswith("old[")]
    ref_names = sorted(set(refs.values()))
    reselect = min_ic >= args.reselect_min_ic
    schemes = [
        ("proteins", ky.bootstrap_weights(len(sizes), args.n_bootstrap, args.seed))
    ]
    if clusters is not None:
        schemes.append(
            (
                "clusters",
                ky.cluster_bootstrap_weights(clusters, args.n_bootstrap, args.seed),
            )
        )
    out = []
    for scheme, weights in schemes:
        pool_old = olds if reselect else ref_names
        wanted = {n: rankings[n] for n in news + pool_old}
        reps = ky.bootstrap(wanted, sizes, weights=weights)
        for name in news + ref_names:
            for metric in BOOT_METRICS:
                ref = refs[metric]
                v, rv = reps[name][metric], reps[ref][metric]
                r_lo, r_hi = ky.ratio_interval(v, rv)
                d_lo, d_hi = np.percentile(v - rv, [2.5, 97.5])
                row = {
                    "aspect": aspect,
                    "min_ic": min_ic,
                    "resampling": scheme,
                    "method": name,
                    "reference": ref,
                    "metric": metric,
                    "ci_lo": float(np.percentile(v, 2.5)),
                    "ci_hi": float(np.percentile(v, 97.5)),
                    "ratio_ci_lo": r_lo,
                    "ratio_ci_hi": r_hi,
                    "diff_ci_lo": float(d_lo),
                    "diff_ci_hi": float(d_hi),
                }
                if reselect:
                    best = np.max([reps[o][metric] for o in olds], axis=0)
                    rr = _ratio(v, best)
                    lo, hi = ky.ratio_interval(v, best)
                    dr_lo, dr_hi = np.percentile(v - best, [2.5, 97.5])
                    row |= {
                        "diff_reselect_ci_lo": float(dr_lo),
                        "diff_reselect_ci_hi": float(dr_hi),
                        "ratio_reselect_ci_lo": lo,
                        "ratio_reselect_ci_hi": hi,
                        "p_ratio_reselect_ge10": float(np.mean(rr >= 10))
                        if rr.size
                        else float("nan"),
                    }
                if metric == "pooled_auprc" and "naive" in reps:
                    lo, hi = ky.ratio_interval(v, reps["naive"][metric])
                    row |= {"ratio_naive_ci_lo": lo, "ratio_naive_ci_hi": hi}
                out.append(row)
        # Precision at the frozen development cutoff, for every method that
        # has one (old variants included, so the deployable comparison has
        # an interval on both sides).
        for name in news + (olds if reselect else ref_names):
            tau = (frozen_tau or {}).get(name)
            if tau is not None and np.isfinite(tau):
                prec = ky.precision_at(rankings[name], tau, weights)
                prec = prec[np.isfinite(prec)]
                if prec.size == 0:  # the cutoff selects nothing in this cohort
                    continue
                out.append(
                    {
                        "aspect": aspect,
                        "min_ic": min_ic,
                        "resampling": scheme,
                        "method": name,
                        "reference": "",
                        "metric": "precision_at_frozen",
                        "ci_lo": float(np.percentile(prec, 2.5)),
                        "ci_hi": float(np.percentile(prec, 97.5)),
                    }
                )
    return out


def evidence_clusters(t1_gaf: Path, aspect: str, proteins) -> np.ndarray:
    """Cluster labels for ``proteins``: connected by shared t1 evidence papers.

    Two cohort proteins are linked when an experimental annotation of each,
    in ``aspect``, cites the same PubMed id; clusters are the connected
    components. A protein citing no paper is its own cluster.
    """
    import gzip

    from src.goa_parser import EXPERIMENTAL_EVIDENCE, has_not_qualifier

    code = {"BP": "P", "MF": "F", "CC": "C"}[aspect]
    index = {p: i for i, p in enumerate(proteins)}
    parent = list(range(len(proteins)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    first_by_paper: dict[str, int] = {}
    with gzip.open(t1_gaf, "rt") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[8] != code or f[1] not in index:
                continue
            if f[6] not in EXPERIMENTAL_EVIDENCE or has_not_qualifier(f[3]):
                continue
            i = index[f[1]]
            for ref in f[5].split("|"):
                if not ref.startswith("PMID:"):
                    continue
                j = first_by_paper.setdefault(ref, i)
                a, b = find(i), find(j)
                if a != b:
                    parent[a] = b
    return np.array([find(i) for i in range(len(proteins))])


def write_outputs(args, rows, boot_rows) -> None:
    from loguru import logger

    args.output_dir.mkdir(parents=True, exist_ok=True)
    metrics = pd.DataFrame(rows)
    first = ["aspect", "min_ic", "method"]
    metrics = metrics[first + [c for c in metrics.columns if c not in first]]
    metrics.to_csv(args.output_dir / "projection_metrics.tsv", sep="\t", index=False)
    if boot_rows:
        pd.DataFrame(boot_rows).to_csv(
            args.output_dir / "projection_bootstrap.tsv", sep="\t", index=False
        )
    (args.output_dir / "projection_args.json").write_text(
        json.dumps({k: str(v) for k, v in vars(args).items()}, indent=2) + "\n"
    )
    logger.info(f"✓ {args.output_dir / 'projection_metrics.tsv'}")


def main(argv=None) -> int:  # pragma: no cover - I/O wiring
    args = parse_args(argv)
    if args.mode == "temporal":
        if args.t1_gaf is None:
            raise SystemExit("--t1-gaf is required in temporal mode")
        return run_temporal(args)
    return run_cv(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
