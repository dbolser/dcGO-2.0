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

    # A named bundle, e.g. the published-dcGO comparison inputs (VALIDATION_PLAN §3)
    uv run python scripts/download_data.py --group dcgo-reference

    # Everything defined in settings.py, including the large optional sources
    uv run python scripts/download_data.py --all

    # List what's available without downloading
    uv run python scripts/download_data.py --list

Notes
-----
* Existing complete files are skipped. Use --force to re-download.
* ``interpro_mappings`` (protein2ipr.dat.gz) is ~20 GB — this can take a while.
* Sources pinned to an immutable release URL carry a ``checksum`` (and
  sometimes a ``size_bytes``) in settings.py — ``disease_ontology``, and the
  frozen published-dcGO and SCOP 1.75 archives. These are verified every time,
  including when an existing file is skipped, so a corrupted, truncated or
  swapped input fails the download step instead of quietly changing a run's
  results.
* A few sources set ``subdir`` and therefore share one directory — the three
  dcGO tables land in ``data/raw/dcgo_reference/``, SCOP in ``data/raw/scop/``.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from urllib.parse import urlparse

import requests

# Make ``config`` importable when run as a plain script from the repo root.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config.settings import ALL_SPECIES_ALIASES, ConfigurationError, Config  # noqa: E402

# The datasets a standard human run actually needs.
DEFAULT_DATASETS = ["goa_annotations", "go_ontology", "interpro_mappings"]

# Named groups, so a reader does not have to know which five sources make up
# "the §3 comparison inputs". Members are ordinary settings.py sources.
DATASET_GROUPS: dict[str, list[str]] = {
    # VALIDATION_PLAN §3: the published dcGO tables plus the SCOP 1.75 release
    # our SSF signatures resolve against. ~135 MB total. Used by
    # validation/compare_original_dcgo.py alongside a --domain-key ssf run.
    "dcgo-reference": [
        "dcgo_domain2go_sql",
        "dcgo_sp2go",
        "dcgo_domain2go_flat",
        "scop_des",
        "scop_hie",
    ],
    # run_dcgo_human.py --ontology hpo needs both files (annotations for the
    # run, hp.obo for --enable-true-path / --enable-relative-inference).
    "hpo": [
        "hpo_annotations",
        "hpo_ontology",
    ],
    # Model-organism phenotype layers: annotations + gene→UniProt id-mapping +
    # the OBO each needs for --enable-true-path / --enable-relative-inference.
    # Run them with --species mouse/worm/zebrafish/fly respectively (see
    # scripts/extract_species_interpro.py for the matching domain universe).
    # Each group carries the annotation files, the gene→UniProt mapping, the
    # OBO, and the per-organism idmapping file that
    # scripts/extract_species_interpro.py builds the species' domain universe
    # from — everything the documented per-species chain needs.
    "mp": [
        "mgi_genepheno",
        "mgi_marker_swissprot",
        "mp_ontology",
        "mouse_idmapping",
    ],
    "wbphenotype": [
        "wormbase_phenotype",
        "worm_idmapping",
        "wbphenotype_ontology",
    ],
    "zfa": [
        "zfin_phenotype",
        "zfin_uniprot",
        "zfa_ontology",
        "zebrafish_idmapping",
    ],
    # fbcv and fbbt share the same three FlyBase tables; only the OBO differs.
    "flybase-phenotype": [
        "flybase_genotype_phenotype",
        "flybase_fbal_to_fbgn",
        "flybase_fbgn_uniprot",
        "fbbt_ontology",
        "fbcv_ontology",
        "fly_idmapping",
    ],
}

# Dated GOA snapshots for the temporal benchmark (VALIDATION_PLAN §2) live in the
# EBI archive, one numbered release per file. The base URL is sourced from
# config/settings.py (single source of truth for dataset locations); we save
# them under a dedicated directory so they never shadow the current-release
# goa_human.gaf.gz.
GOA_ARCHIVE_DIRNAME = "goa_archive"

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


class ChecksumError(RuntimeError):
    """A downloaded file did not match the checksum pinned in settings.py."""


def file_digest(path: Path, algorithm: str) -> str:
    """Hex digest of ``path``, streamed so large files never load into memory."""
    digest = hashlib.new(algorithm)
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_checksum(path: Path, expected: tuple[str, str]) -> None:
    """Check ``path`` against an ``(algorithm, hex digest)`` pair.

    Only pinned (immutable-URL) sources carry a checksum, and this is what makes
    the pinning worth anything: a run's inputs are then reproducible by content,
    not merely by URL.

    Raises:
        ChecksumError: on a mismatch.
    """
    algorithm, want = expected
    got = file_digest(path, algorithm)
    if got != want:
        raise ChecksumError(
            f"{algorithm} mismatch for {path}: expected {want}, got {got}"
        )
    print(f"  ✔ {algorithm} verified: {want}")


def verify_size(path: Path, size_bytes: int) -> None:
    """Check ``path`` against the byte count pinned in settings.py.

    Cheaper than a digest and it catches the common failure directly: a
    truncated download. Carried by the frozen archives (the published dcGO
    tables, SCOP 1.75), where a short read would silently change a §3
    comparison rather than fail it.

    Raises:
        ChecksumError: on a mismatch, so callers handle one exception type.
    """
    actual = path.stat().st_size
    if actual != size_bytes:
        raise ChecksumError(
            f"{path.name}: expected {size_bytes:,} bytes, got {actual:,}. "
            "Truncated download, or the source changed."
        )
    print(f"  ✔ size verified: {size_bytes:,} bytes")


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

    # Clean up the partial file on any interruption (network drop, Ctrl-C, disk
    # full) so a failed ~20 GB download doesn't silently linger. We re-raise so
    # the caller still sees the failure.
    try:
        with requests.get(url, stream=True, timeout=timeout) as resp:
            resp.raise_for_status()
            # Content-Length may be absent or malformed; treat as "unknown size".
            try:
                total = int(resp.headers.get("Content-Length") or 0)
            except ValueError:
                total = 0
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
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise

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
        "--group",
        action="append",
        choices=sorted(DATASET_GROUPS),
        metavar="NAME",
        help="Download a named bundle of datasets (repeatable). Available: "
        + ", ".join(sorted(DATASET_GROUPS)),
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
        "--species",
        default="human",
        metavar="NAME",
        help="Species whose GOA annotations to download (default: human). Any "
        "other value (e.g. 'mouse', 'zebrafish') retargets the goa_annotations "
        "dataset to <base>/<SPECIES_UPPER>/goa_<species>.gaf.gz and saves it as "
        "goa_annotations/goa_<species>.gaf.gz.",
    )
    parser.add_argument(
        "--goa-archive",
        action="append",
        metavar="VERSION",
        help="Download a dated human GOA snapshot (numbered EBI release, e.g. "
        "205 = 2021-04) into data/raw/goa_archive/ for the temporal benchmark "
        "(VALIDATION_PLAN §2). Repeatable.",
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
        print("\nGroups (--group NAME):\n")
        for group, members in sorted(DATASET_GROUPS.items()):
            print(f"  {group:22s} {', '.join(members)}")
        return 0

    raw_dir = config.DATA_DIR / "raw"

    # Dated GOA snapshots are handled separately from the settings.py sources:
    # they share one archive URL pattern and land in their own directory.
    if args.goa_archive:
        print(f"Saving dated GOA snapshots under: {raw_dir / GOA_ARCHIVE_DIRNAME}\n")
        archive_failures: list[str] = []
        for version in args.goa_archive:
            filename = f"goa_human.gaf.{version}.gz"
            url = f"{config.goa_archive_base_url}/{filename}"
            dest = raw_dir / GOA_ARCHIVE_DIRNAME / filename
            print(f"[goa_archive {version}] dated human GOA snapshot")
            try:
                download_one(url, dest, force=args.force, timeout=args.timeout)
            except (requests.RequestException, OSError) as exc:
                print(f"  ✘ FAILED: {exc}")
                archive_failures.append(version)
            print()
        if archive_failures:
            print(f"GOA archive download errors. Failed: {', '.join(archive_failures)}")
            return 1
        if not (args.all or args.datasets or args.group):
            # Only archive snapshots were requested — done.
            print("All requested GOA snapshots are in place.")
            return 0

    if args.all:
        selected = list(sources.keys())
    elif args.datasets or args.group:
        # De-duplicate while preserving order, so --group dcgo-reference
        # --datasets scop_hie does not fetch scop_hie twice.
        selected = []
        for name in list(args.datasets or []) + [
            member for group in (args.group or []) for member in DATASET_GROUPS[group]
        ]:
            if name not in selected:
                selected.append(name)
    else:
        selected = DEFAULT_DATASETS

    print(f"Saving datasets under: {raw_dir}\n")

    species = (args.species or "human").strip().lower()

    failures: list[str] = []
    for name in selected:
        ds = sources[name]
        url = ds.url
        dest = raw_dir / (ds.subdir or name) / _filename_for(ds.url, name)
        description = ds.description

        # GOA is per-species; retarget the pinned (human) URL for any other
        # organism so a fresh checkout can go download → extract → run for e.g.
        # mouse without editing settings.py.
        if name == "goa_annotations" and species != "human":
            try:
                url = config.goa_url_for(species)
            except ConfigurationError as exc:
                # Temporal snapshots have no URL composable from the name; say
                # so and point at the flag that does fetch them, rather than
                # tracebacking or downloading a 404.
                print(f"[{name}] ✘ {exc}")
                failures.append(name)
                print()
                continue
            # Name the destination after the REQUESTED species, not after the
            # URL. For every organism the two agree (MOUSE/goa_mouse.gaf.gz),
            # but the all-species background is served from one cross-organism
            # file — goa_uniprot_all.gaf.gz — while the pipeline looks for
            # goa_allspecies.gaf.gz. Deriving the filename from the URL there
            # would download the right bytes to a path no later step reads.
            dest = raw_dir / name / f"goa_{species}.gaf.gz"
            description = f"Gene Ontology Annotation (GOA) database for {species}"
            if species in ALL_SPECIES_ALIASES:
                print(
                    f"[{name}] note: this is the full cross-organism release "
                    "(~11.7 GB, ~97% IEA). run_dcgo_human.py applies its own "
                    "evidence filter, so it is usable as-is, but every run then "
                    "re-parses all of it. scratch_allspecies/01_build_gaf.sh "
                    "pre-filters it to non-IEA once (~258 MB)."
                )

        print(f"[{name}] {description}")
        try:
            download_one(url, dest, force=args.force, timeout=args.timeout)
            # Verify whether or not we just fetched it: these are only set
            # for immutable-release URLs, and an already-present file is exactly
            # the case where silent corruption would go unnoticed.
            if (size_bytes := getattr(ds, "size_bytes", None)) is not None:
                verify_size(dest, size_bytes)
            if (expected := ds.checksum_parts()) is not None:
                verify_checksum(dest, expected)
        except (requests.RequestException, OSError, ChecksumError) as exc:
            print(f"  ✘ FAILED: {exc}")
            # Release-stamped URLs (WormBase, FlyBase) break on upstream
            # release turnover; their DataSource says what to update.
            if ds.update_hint:
                print(f"    hint: {ds.update_hint}")
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
