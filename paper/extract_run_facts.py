"""Copy the counts a run only logs into a result file next to its metrics.

The benchmark driver logs, but does not tabulate, a few counts the paper
cites: background size, same-gene accessions removed, cohort proteins with
no t0 record, sequence and literature coverage, and cohort sizes. This
reads them from the run's log and writes ``<run>/run_facts.tsv``, keeping
the log line each value came from.

    uv run python paper/extract_run_facts.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RUNS = {
    "dev2018_frozen": ROOT / "logs/projection/dev2018_frozen.log",
    "final": ROOT / "logs/projection/final.log",
}
#: The literature features each run read, and the t0 GAF whose coverage counts.
LITERATURE = {
    "dev2018_frozen": ("medembed_protein_y2017.json", "goa_human.gaf.181.gz"),
    "final": ("medembed_protein_y2020.json", "goa_human.gaf.205.gz"),
}
PATTERNS = {
    "background_annotated_proteins": r"Background: ([\d,]+) annotated proteins",
    "same_gene_accessions_dropped": r"\(([\d,]+) same-gene human accessions dropped\)",
    "background_proteins_with_cohort_unit": r"dropped\), ([\d,]+) with a cohort unit",
    "cohort_proteins_without_t0_record": r"no t0 record \(no embedding used\): ([\d,]+)",
    "sequence_embedding_proteins": r"Sequence embeddings: ([\d,]+) proteins",
    "literature_vector_accessions": r"literature ([\d,]+) proteins",
    "cohort_BP": r"\[BP\] ([\d,]+) cohort proteins",
    "cohort_MF": r"\[MF\] ([\d,]+) cohort proteins",
    "cohort_CC": r"\[CC\] ([\d,]+) cohort proteins",
}


def main() -> int:
    for run, log in RUNS.items():
        text = log.read_text(encoding="utf-8")
        rows = []
        for key, pattern in PATTERNS.items():
            m = re.search(pattern, text)
            if not m:
                raise SystemExit(f"{run}: {key} not found in {log}")
            line = text[text.rfind("\n", 0, m.start()) + 1 : text.find("\n", m.end())]
            rows.append(
                {
                    "fact": key,
                    "value": int(m.group(1).replace(",", "")),
                    "log_line": line.split(" - ", 1)[-1].strip(),
                }
            )
        sidecar, gaf = LITERATURE[run]
        meta = json.loads(
            (ROOT / "data/processed/literature" / sidecar).read_text(encoding="utf-8")
        )
        cover = meta["coverage"][gaf]
        for key, value, source in (
            ("literature_papers_up_to_cutoff", meta["counts"]["papers_le_cutoff"],
             f"{sidecar}: counts.papers_le_cutoff"),
            ("t0_gaf_accessions", cover["accessions"], f"{sidecar}: coverage.{gaf}"),
            ("t0_gaf_accessions_with_papers", cover["with_papers"],
             f"{sidecar}: coverage.{gaf}"),
        ):  # fmt: skip
            rows.append({"fact": key, "value": int(value), "log_line": source})
        out = ROOT / "validation/projection" / run / "run_facts.tsv"
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        print(out)
        print(pd.DataFrame(rows)[["fact", "value"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
