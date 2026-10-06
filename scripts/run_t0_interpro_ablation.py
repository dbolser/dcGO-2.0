"""Rerun the §4 ``ipr_manual`` ablation cell on InterPro 85.0 (t0) architectures.

TODO.md P0 "Remove temporal look-ahead from domain architectures": the committed
cell trains on GOA release 205 (2021-04-21) but reads the 2026-07 protein2ipr,
both in training and in the evaluator's transfer step. InterPro 85.0
(2021-04-08) is the last release before the t0 GAF. Three steps, each printed
before it runs and logged to ``<run-dir>/<step>.log``:

1. **Extract** the human subset of the archived protein2ipr. Proteins are
   selected by the t0 *and* t1 GAFs under any evidence, so every training
   protein and every held-out protein the evaluator transfers to has t0
   architectures. (The current subset was selected by the t1 GAF's non-IEA
   proteins alone, a strict subset of this.) Skipped when the subset exists
   and its provenance marker names the same source.
2. **Train** the nine rungs: exactly the commands recorded in
   ``validation/ablation_manifests/<rung>.json``, with ``--output-dir`` moved
   under ``--run-dir`` and ``--gaf``/``--interpro`` added.
3. **Evaluate** with ``validation/ablation.py`` under the committed
   evaluation's arguments (all defaults besides the GAFs and run dir), plus
   ``--interpro`` set to the subset. The rung manifests are copied to
   ``<eval-dir>/manifests/``, as for the other cells.

Run it from the repository root; paths are relative to it, and ``data/`` must
hold the inputs. The md5 of the archived protein2ipr is verified first
(``scripts/download_data.py --group interpro-85``). Resumable: a rung whose
manifest says ``completed`` is not rerun — safe, because the evaluator refuses
a rung trained on a different protein2ipr than ``--subset``.

Dry run of the whole chain without the archive: pass the current subset as
``--source`` (with scratch ``--subset``/``--run-dir``/``--eval-dir``). A wider
selection filtered out of a file that holds only selected proteins is that
file again, so it must reproduce the committed ``ipr_manual`` numbers exactly.

    uv run python scripts/run_t0_interpro_ablation.py
    uv run python scripts/run_t0_interpro_ablation.py --dry-run   # print only
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.run_manifest import manifest_filename  # noqa: E402
from src.universe_provenance import read_marker  # noqa: E402

MANIFESTS = REPO / "validation" / "ablation_manifests"
ARCHIVE = Path("data/raw/interpro_archive/85.0/protein2ipr.dat.gz")
T0_GAF = Path("data/raw/goa_archive/goa_human.gaf.205.gz")
T1_GAF = Path("data/raw/goa_annotations/goa_human.gaf.gz")
RUNGS = (
    "single",
    "supra",
    "supra_input",
    "supra_relative",
    "supra_output",
    "supra_input_relative",
    "supra_input_output",
    "supra_relative_output",
    "full",
)


def rung_command(manifest: Path, output_dir: Path, interpro: Path) -> list[str]:
    """The command a rung's manifest records, re-pointed at t0 inputs.

    Only the interpreter (``command[0]``), the ``--output-dir`` value and the
    two input flags change; every analysis flag is the recorded one.
    """
    args = json.loads(manifest.read_text(encoding="utf-8"))["command"][1:]
    args[args.index("--output-dir") + 1] = str(output_dir)
    return [sys.executable, *args, "--gaf", str(T0_GAF), "--interpro", str(interpro)]


def completed(output_dir: Path) -> bool:
    """Whether a rung's run finished (its manifest is finalised)."""
    manifest = output_dir / manifest_filename("go")
    return (
        manifest.exists()
        and json.loads(manifest.read_text(encoding="utf-8"))["status"] == "completed"
    )


def build_steps(
    source: Path, subset: Path, run_dir: Path, eval_dir: Path
) -> list[tuple[str, list[str]]]:
    """(log name, command) for every step still to run, in order."""
    steps = []
    if source == ARCHIVE:
        steps.append(
            (
                "verify",
                [sys.executable, "scripts/download_data.py", "--group", "interpro-85"],
            )
        )
    marker = read_marker(subset)
    if not (subset.exists() and marker and marker.interpro_source == str(source)):
        steps.append(
            (
                "extract",
                [
                    sys.executable,
                    "extract_human_interpro.py",
                    "--source",
                    str(source),
                    "--gaf",
                    str(T0_GAF),
                    "--gaf",
                    str(T1_GAF),
                    "--evidence-filter",
                    "all",
                    "--output",
                    str(subset),
                ],
            )
        )
    for rung in RUNGS:
        if not completed(run_dir / rung):
            steps.append(
                (rung, rung_command(MANIFESTS / f"{rung}.json", run_dir / rung, subset))
            )
    steps.append(
        (
            "eval",
            [
                sys.executable,
                "validation/ablation.py",
                "--t0-gaf",
                str(T0_GAF),
                "--t1-gaf",
                str(T1_GAF),
                "--run-dir",
                str(run_dir),
                "--interpro",
                str(subset),
                "--output-dir",
                str(eval_dir),
            ],
        )
    )
    return steps


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Rerun the ipr_manual ablation cell on InterPro 85.0 (t0) "
        "domain architectures (VALIDATION_PLAN §4)."
    )
    parser.add_argument("--source", type=Path, default=ARCHIVE)
    parser.add_argument(
        "--subset",
        type=Path,
        default=Path("data/interim/protein2ipr_human_t0_ipr85.dat.gz"),
    )
    parser.add_argument(
        "--run-dir", type=Path, default=Path("results/ablation-t0interpro/ipr_manual")
    )
    parser.add_argument(
        "--eval-dir",
        type=Path,
        default=Path("validation/ablation_cells/ipr_manual_t0interpro"),
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="print the commands and exit"
    )
    args = parser.parse_args(argv)

    steps = build_steps(args.source, args.subset, args.run_dir, args.eval_dir)
    if args.dry_run:
        for name, cmd in steps:
            print(f"# {name}\n{' '.join(cmd)}")
        return 0

    args.run_dir.mkdir(parents=True, exist_ok=True)
    for name, cmd in steps:
        log = args.run_dir / f"{name}.log"
        print(f"[{name}] {' '.join(cmd)}  > {log}", flush=True)
        with open(log, "w") as handle:
            result = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            print(f"[{name}] FAILED rc={result.returncode} — see {log}", flush=True)
            return 1

    manifests = args.eval_dir / "manifests"
    manifests.mkdir(parents=True, exist_ok=True)
    for rung in RUNGS:
        shutil.copy(
            args.run_dir / rung / "run_manifest_go.json", manifests / f"{rung}.json"
        )
    print(f"Done. Evaluation tables and rung manifests in {args.eval_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
