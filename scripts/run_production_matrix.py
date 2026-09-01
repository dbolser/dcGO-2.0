"""Run the production matrix: every registry ontology, single + supra domains.

The deliverable driver for the 2026-08 expansion (see data/ACQUISITION_MATRIX.md
and the project mandate): one run per ontology layer in each applicable
configuration, sequentially, with manifests, restartable.

Configurations per layer:

* ``baseline`` — the overall inference alone, every layer. For flat
  cross-reference vocabularies this is the only possible configuration.
* ``paperparity`` — hierarchical layers only: input True-Path propagation,
  relative inference combined before BH, True-Path propagation of the results,
  and the ``--min-ic 1`` reporting floor (the recommended production floor;
  kills the DAG roots and the near-universal band, VALIDATION_PLAN item 2).

GO additionally runs with ``--species allspecies`` (the winning training
universe, MULTISPECIES_BACKGROUND.md) and in ``--evidence-filter experimental``
mode (GAF evidence codes; the direct-vs-IEA contrast from the mandate). Model
organism layers run under their own ``--species`` — the dcGO trick: the
association is learned on that organism's proteins, domains are
species-agnostic.

Restartable: a run whose output directory already contains its run manifest is
skipped, so a crashed or interrupted matrix resumes where it stopped. Per-run
failures are recorded and do not stop the matrix.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.ontology_registry import ONTOLOGIES  # noqa: E402
from src.run_manifest import manifest_filename  # noqa: E402

#: Layers learned on a non-human organism's own proteins.
SPECIES_FOR_ONTOLOGY = {
    "mp": "mouse",
    "wbphenotype": "worm",
    "wbbt": "worm",
    "zfa": "zebrafish",
    "fbcv": "fly",
    "fbbt": "fly",
}

#: Registry keys that are not standalone production layers.
EXCLUDED = {
    "xref",  # generic escape hatch; needs --xref-db, not a curated layer
}

PAPER_PARITY_FLAGS = [
    "--propagate-annotations",
    "--enable-relative-inference",
    "--enable-true-path",
    "--min-ic",
    "1",
]


def build_matrix() -> list[tuple[str, list[str]]]:
    """(run_name, extra CLI args) for every cell of the matrix."""
    cells: list[tuple[str, list[str]]] = []
    for key, entry in sorted(ONTOLOGIES.items()):
        if key in EXCLUDED:
            continue
        species = SPECIES_FOR_ONTOLOGY.get(key)
        base = ["--ontology", key]
        name = key
        if species is not None:
            base += ["--species", species]
            name = f"{key}_{species}"
        cells.append((f"{name}_baseline", list(base)))
        if entry.supports_true_path:
            cells.append((f"{name}_paperparity", base + PAPER_PARITY_FLAGS))

    # GO extras: the allspecies training universe, and the experimental-evidence
    # (direct, non-IEA) contrast, in both configurations.
    go = ["--ontology", "go"]
    cells.append(("go_allspecies_baseline", go + ["--species", "allspecies"]))
    cells.append(
        (
            "go_allspecies_paperparity",
            go + ["--species", "allspecies"] + PAPER_PARITY_FLAGS,
        )
    )
    cells.append(
        ("go_experimental_baseline", go + ["--evidence-filter", "experimental"])
    )
    cells.append(
        (
            "go_experimental_paperparity",
            go + ["--evidence-filter", "experimental"] + PAPER_PARITY_FLAGS,
        )
    )
    return cells


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=Path("results/production"))
    parser.add_argument("--num-cores", type=int, default=8)
    parser.add_argument(
        "--only", nargs="*", default=None, help="run only these cell names"
    )
    parser.add_argument("--list", action="store_true", help="print the matrix and exit")
    args = parser.parse_args(argv)

    cells = build_matrix()
    if args.only:
        wanted = set(args.only)
        unknown = wanted - {name for name, _ in cells}
        if unknown:
            parser.error(f"unknown cell names: {sorted(unknown)}")
        cells = [(name, extra) for name, extra in cells if name in wanted]
    if args.list:
        for name, extra in cells:
            print(f"{name}\t{' '.join(extra)}")
        return 0

    args.output_root.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []
    skipped = 0
    for i, (name, extra) in enumerate(cells, 1):
        out_dir = args.output_root / name
        ontology = extra[extra.index("--ontology") + 1]
        if (out_dir / manifest_filename(ontology)).exists():
            print(f"[{i}/{len(cells)}] {name}: manifest exists, skipping")
            skipped += 1
            continue
        cmd = [
            sys.executable,
            "run_dcgo_human.py",
            *extra,
            "--num-cores",
            str(args.num_cores),
            "--output-dir",
            str(out_dir),
        ]
        log_path = args.output_root / f"{name}.log"
        print(f"[{i}/{len(cells)}] {name}: running", flush=True)
        start = time.time()
        with open(log_path, "w") as log:
            result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        wall = time.time() - start
        if result.returncode != 0:
            failures.append(name)
            print(
                f"[{i}/{len(cells)}] {name}: FAILED rc={result.returncode} "
                f"({wall:.0f}s) — see {log_path}",
                flush=True,
            )
        else:
            print(f"[{i}/{len(cells)}] {name}: done ({wall:.0f}s)", flush=True)

    print(
        f"\nMatrix complete: {len(cells) - len(failures) - skipped} ran, "
        f"{skipped} skipped (already done), {len(failures)} failed"
    )
    if failures:
        print("Failed cells: " + ", ".join(failures))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
