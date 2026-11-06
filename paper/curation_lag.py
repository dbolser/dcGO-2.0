"""How much of the benchmark truth cites papers published before t0?

For each split and GO aspect, reads the t1 GAF's experimental lines for the
no-knowledge cohort (as built by ``validation/ablation.py``), takes the
PubMed ids they cite, and reports the share of evidence lines and of
distinct papers published up to the literature cutoff year (the year before
t0). A high share means the later experimental annotation rests on papers
that already existed at t0 (curation lag), which the literature channel can
see. Writes ``validation/projection/<run>/curation_lag.tsv``.

    uv run --with pyarrow python paper/curation_lag.py
"""

from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
_spec = importlib.util.spec_from_file_location(
    "ablation", ROOT / "validation" / "ablation.py"
)
ab = importlib.util.module_from_spec(_spec)
sys.modules["ablation"] = ab
_spec.loader.exec_module(ab)

from src.goa_parser import EXPERIMENTAL_EVIDENCE, has_not_qualifier  # noqa: E402

ABSTRACTS = Path("/home/danbolser/Work/BioEmbed/genesearch/index/abstracts.parquet")
SPLITS = {
    "dev2018_frozen": dict(
        t0="data/raw/goa_archive/goa_human.gaf.181.gz",
        t1="data/raw/goa_archive/goa_human.gaf.205.gz",
        interpro="data/interim/protein2ipr_human_t0_2018_ipr69.dat.gz",
        cohort="data/interim/protein2ipr_human_t0_ipr85.dat.gz",
        cutoff=2017,
    ),
    "final": dict(
        t0="data/raw/goa_archive/goa_human.gaf.205.gz",
        t1="data/raw/goa_annotations/goa_human.gaf.gz",
        interpro="data/interim/protein2ipr_human_t0_ipr85.dat.gz",
        cohort="data/interim/protein2ipr_human.dat.gz",
        cutoff=2020,
    ),
}
ASPECT = {"P": "BP", "F": "MF", "C": "CC"}


def main() -> int:
    years = (
        pq.read_table(ABSTRACTS, columns=["pmid", "year"])
        .to_pandas()
        .set_index("pmid")["year"]
    )
    for run, s in SPLITS.items():
        inputs = ab.load_benchmark_inputs(
            ROOT / s["t0"],
            ROOT / s["t1"],
            ROOT / s["interpro"],
            ROOT / "data/raw/go_ontology/go-basic.obo",
            "interpro",
            "manual",
            cohort_interpro=ROOT / s["cohort"],
        )

        # Two views: every experimental line of the cohort, and only the lines
        # behind the scored truth (proteins and terms at IC >= 2: a line
        # counts if its term or one of its ancestors has IC >= 2).
        def informative(term: str) -> bool:
            term = inputs.alt_id_map.get(term, term)
            term = inputs.replaced_by_map.get(term, term)
            closure = {term} | inputs.get_ancestors(term)
            return any(inputs.ic.get(t, 0.0) >= 2.0 for t in closure)

        scored = {
            a: {p for p, terms in inputs.benchmark[a].items()
                if any(inputs.ic.get(t, 0.0) >= 2.0 for t in terms)}
            for a in ASPECT.values()
        }  # fmt: skip
        views = {"all_terms": {}, "ic2": {}}
        for v in views.values():
            v.update({a: [] for a in ASPECT.values()})
        with gzip.open(ROOT / s["t1"], "rt") as handle:
            for line in handle:
                if line.startswith("!"):
                    continue
                f = line.split("\t")
                if f[6] not in EXPERIMENTAL_EVIDENCE or has_not_qualifier(f[3]):
                    continue
                aspect = ASPECT.get(f[8])
                if not aspect or f[1] not in inputs.benchmark[aspect]:
                    continue
                pmids = [
                    int(r[5:])
                    for r in f[5].split("|")
                    if r.startswith("PMID:") and r[5:].isdigit()
                ]
                views["all_terms"][aspect].extend(pmids)
                if f[1] in scored[aspect] and informative(f[4]):
                    views["ic2"][aspect].extend(pmids)
        rows = []
        for (view, aspect), pmids in (
            ((view, a), pm) for view, by in views.items() for a, pm in by.items()
        ):
            arr = np.array(pmids)
            y = years.reindex(arr).to_numpy()
            distinct = np.unique(arr)
            yd = years.reindex(distinct).to_numpy()
            rows.append(
                {
                    "view": view,
                    "aspect": aspect,
                    "cutoff_year": s["cutoff"],
                    "evidence_lines": len(arr),
                    "share_lines_pre_cutoff": float(np.nanmean(y <= s["cutoff"]))
                    if len(arr)
                    else float("nan"),
                    "distinct_papers": len(distinct),
                    "share_papers_pre_cutoff": float(np.nanmean(yd <= s["cutoff"]))
                    if len(distinct)
                    else float("nan"),
                    "lines_without_year": int(np.isnan(y.astype(float)).sum()),
                }
            )
        out = ROOT / "validation" / "projection" / run / "curation_lag.tsv"
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        print(run)
        print(pd.DataFrame(rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
