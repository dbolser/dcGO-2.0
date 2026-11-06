"""Replay a committed ablation cell's rung commands on other inputs.

Each ``<rung>.json`` run manifest in ``--manifests`` records the exact
``run_dcgo_human.py`` command of one §4 rung. This reruns those commands with
only the inputs and output directory changed — ``--gaf``, ``--interpro`` and
``--output-dir`` — so a different split (e.g. the 2018 development split of
``KNOWLEDGE_PROJECTION.md``) is trained with the same analysis flags.
A rung whose manifest in ``--run-dir`` already says ``completed`` with the
same input hashes is skipped.

    uv run python scripts/replay_ladder.py \\
        --manifests validation/ablation_cells/ipr_manual_t0interpro_v2/manifests \\
        --gaf data/raw/goa_archive/goa_human.gaf.181.gz \\
        --interpro data/interim/protein2ipr_human_t0_2018_ipr69.dat.gz \\
        --run-dir results/projection-dev2018/ipr_manual
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.run_manifest import manifest_filename, sha256_file  # noqa: E402

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


def replay_command(manifest: Path, gaf: Path, interpro: Path, out: Path) -> list[str]:
    """The recorded command with its inputs and output dir replaced."""
    args = json.loads(manifest.read_text(encoding="utf-8"))["command"][1:]
    for flag, value in (
        ("--output-dir", out),
        ("--gaf", gaf),
        ("--interpro", interpro),
    ):
        if flag in args:
            args[args.index(flag) + 1] = str(value)
        else:
            args += [flag, str(value)]
    return [sys.executable, *args]


def up_to_date(out: Path, gaf_sha: str, interpro_sha: str) -> bool:
    path = out / manifest_filename("go")
    if not path.exists():
        return False
    m = json.loads(path.read_text(encoding="utf-8"))
    hashes = json.dumps(m.get("inputs", m))
    return (
        m.get("status") == "completed" and gaf_sha in hashes and interpro_sha in hashes
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--manifests", type=Path, required=True)
    parser.add_argument("--gaf", type=Path, required=True)
    parser.add_argument("--interpro", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--rungs", nargs="+", default=list(RUNGS))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    gaf_sha, ipr_sha = sha256_file(args.gaf), sha256_file(args.interpro)
    for rung in args.rungs:
        out = args.run_dir / rung
        if up_to_date(out, gaf_sha, ipr_sha):
            print(f"[{rung}] up to date, skipped")
            continue
        cmd = replay_command(
            args.manifests / f"{rung}.json", args.gaf, args.interpro, out
        )
        print(f"[{rung}] {' '.join(cmd[1:])}", flush=True)
        if args.dry_run:
            continue
        out.mkdir(parents=True, exist_ok=True)
        with open(out.parent / f"{rung}.log", "w", encoding="utf-8") as log:
            result = subprocess.run(cmd, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            print(f"[{rung}] FAILED (exit {result.returncode}); see {rung}.log")
            return result.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
