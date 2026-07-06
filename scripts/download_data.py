#!/usr/bin/env python
"""Download the raw datasets required for a real dcGO run.

This is a simple, dependency-light downloader (uses ``requests``). It reads the
dataset URLs from ``config/settings.py`` so there is a single source of truth,
and saves each file to the layout the rest of the pipeline expects:

    data/raw/<source_name>/<original_filename>

e.g.

    data/raw/goa_annotations/goa_human.gaf.gz
    data/raw/go_ontology/go-basic.obo
    data/raw/interpro_mappings/protein2ipr.dat.gz

Those are exactly the paths ``extract_human_interpro.py`` and
``run_dcgo_human.py`` look for, so a fresh checkout can go straight from
download → extract → run.

Usage
-----
    # The three required datasets (GOA human, GO ontology, InterPro mappings)
    uv run python scripts/download_data.py

    # A specific subset (repeat --datasets)
    uv run python scripts/download_data.py --datasets goa_annotations --datasets go_ontology

    # Everything defined in settings.py, including the large optional sources
    uv run python scripts/download_data.py --all

    # List what's available without downloading
    uv run python scripts/download_data.py --list

Notes
-----
* Existing complete files are skipped. Use --force to re-download.
* ``interpro_mappings`` (protein2ipr.dat.gz) is ~20 GB — this can take a while.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.parse import urlparse

import requests

# Make ``config`` importable when run as a plain script from the repo root.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config.settings import Config  # noqa: E402

# The datasets a standard human run actually needs.
DEFAULT_DATASETS = ["goa_annotations", "go_ontology", "interpro_mappings"]

CHUNK_SIZE = 1 << 20  # 1 MiB


def _filename_for(url: str, source_name: str) -> str:
    """Original filename from the URL, falling back to the source name."""
    name = Path(urlparse(url).path).name
    return name or f"{source_name}.data"


def _human_size(num_bytes: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if num_bytes < 1024 or unit == "TB":
            return f"{num_bytes:.1f} {unit}"
        num_bytes /= 1024
    return f"{num_bytes:.1f} TB"


def download_one(url: str, dest: Path, *, force: bool, timeout: int) -> bool:
    """Download ``url`` to ``dest``. Returns True if a download happened."""
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and not force:
        print(
            f"  ✔ already present, skipping: {dest} ({_human_size(dest.stat().st_size)})"
        )
        return False

    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"  ↓ {url}\n    → {dest}")

    with requests.get(url, stream=True, timeout=timeout) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        with open(tmp, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=CHUNK_SIZE):
                if not chunk:
                    continue
                fh.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded / total * 100
                    print(
                        f"\r    {_human_size(downloaded)} / {_human_size(total)} ({pct:5.1f}%)",
                        end="",
                        flush=True,
                    )
                else:
                    print(f"\r    {_human_size(downloaded)}", end="", flush=True)
        print()

    tmp.replace(dest)
    print(f"  ✔ done: {dest} ({_human_size(dest.stat().st_size)})")
    return True


def main() -> int:
    config = Config()
    sources = config.data_sources

    parser = argparse.ArgumentParser(
        description="Download raw datasets for the dcGO pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--datasets",
        action="append",
        choices=sorted(sources.keys()),
        metavar="NAME",
        help="Dataset to download (repeatable). Defaults to the required set: "
        + ", ".join(DEFAULT_DATASETS),
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Download every dataset defined in settings.py (includes large optional sources).",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available datasets and exit.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if the file already exists.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="Per-request connect/read timeout in seconds (default: 60).",
    )
    args = parser.parse_args()

    if args.list:
        print("Available datasets:\n")
        for name, ds in sorted(sources.items()):
            flag = "required" if ds.required else "optional"
            print(f"  {name:22s} [{flag}]  {ds.description}")
            print(f"  {'':22s}          {ds.url}")
        return 0

    if args.all:
        selected = list(sources.keys())
    elif args.datasets:
        selected = args.datasets
    else:
        selected = DEFAULT_DATASETS

    raw_dir = config.DATA_DIR / "raw"
    print(f"Saving datasets under: {raw_dir}\n")

    failures: list[str] = []
    for name in selected:
        ds = sources[name]
        dest = raw_dir / name / _filename_for(ds.url, name)
        print(f"[{name}] {ds.description}")
        try:
            download_one(ds.url, dest, force=args.force, timeout=args.timeout)
        except (requests.RequestException, OSError) as exc:
            print(f"  ✘ FAILED: {exc}")
            failures.append(name)
        print()

    if failures:
        print(f"Completed with errors. Failed: {', '.join(failures)}")
        return 1

    print("All requested datasets are in place.")
    print("Next: uv run python extract_human_interpro.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
