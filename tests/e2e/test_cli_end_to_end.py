"""End-to-end tests: run the CLI on fixtures and check the output is usable.

`ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` (P0) asks for "a small end-to-end test
that invokes the installed CLI on fixtures and verifies schema-valid output",
because "the unit and integration tests do not currently demonstrate that the
installed command works". CI already smoke-tests `dcgo --help` against a built
wheel, which proves the entry point resolves but nothing about the pipeline.

These tests drive the real command in a subprocess against a synthetic dataset
small enough to run in a second, and then assert on the output file rather than
on the exit code alone. The fixture is built so that the *correct* answer is
known in advance: `IPR_HIT` is carried by exactly the proteins annotated with
`GO:0000001` and by no others, so a working pipeline must report that pair,
must not report the deliberately uninformative `IPR_NOISE`, and must place the
true pair at the top of the ranking.

The dataset is laid out under a temporary working directory in the layout the
runner expects (`data/interim/protein2ipr_<species>.dat.gz` and
`data/raw/goa_annotations/goa_<species>.gaf.gz`), so this also pins those
conventions: if they change, this test fails rather than a user discovering it.
"""

from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]

SPECIES = "fixture"

#: Columns `run_dcgo_human.py` writes for the significant-association table.
#: A consumer that reads this file by name depends on every one of them.
EXPECTED_COLUMNS = [
    "domain",
    "go_term",
    "p_value",
    "adj_p_value",
    "odds_ratio",
    "hyper_score",
    "domain_type",
    "constituent_domains",
    "n_observations",
]

#: Carried by exactly the GO:0000001-annotated proteins — the pair that must be found.
IPR_HIT = "IPR000001"
#: Carried by exactly half of each group, so its rate matches the background.
IPR_NOISE = "IPR000002"
#: Carried by everyone. Its only job is to put every protein in the universe:
#: a protein with no domain at all is dropped by `restrict_to_universe`, and if
#: the dropped proteins are not label-balanced then `IPR_NOISE` stops being
#: uninformative *within the surviving universe* and is correctly reported.
IPR_UBIQUITOUS = "IPR000003"
TRUE_TERM = "GO:0000001"
OTHER_TERM = "GO:0000002"

N_POSITIVE = 60
N_NEGATIVE = 60


def _gaf_rows(n_positive: int, n_negative: int) -> list[str]:
    """GAF 2.2 rows: positives get GO:0000001, negatives get GO:0000002."""
    rows = ["!gaf-version: 2.2\n"]
    for index in range(n_positive + n_negative):
        accession = f"P{index:05d}"
        term = TRUE_TERM if index < n_positive else OTHER_TERM
        fields = [
            "UniProtKB",
            accession,
            accession,
            "",
            term,
            "PMID:1",
            "IDA",
            "",
            "F",
            "",
            "",
            "protein",
            "taxon:9606",
            "20200101",
            "UniProt",
            "",
            "",
        ]
        rows.append("\t".join(fields) + "\n")
    return rows


def _protein2ipr_rows(n_positive: int, n_negative: int) -> list[str]:
    """Domain rows, ordered along the sequence so supra-domains are well defined.

    `IPR_HIT` tracks the label exactly; `IPR_NOISE` sits on every other protein
    in *both* groups so its rate equals the background; `IPR_UBIQUITOUS` is on
    everything so no protein leaves the universe.
    """
    rows = []
    for index in range(n_positive + n_negative):
        accession = f"P{index:05d}"
        if index < n_positive:
            rows.append(f"{accession}\t{IPR_HIT}\tHit domain\tPF00001\t10\t120\n")
        if index % 2 == 0:
            rows.append(f"{accession}\t{IPR_NOISE}\tNoise domain\tPF00002\t200\t320\n")
        rows.append(f"{accession}\t{IPR_UBIQUITOUS}\tEverywhere\tPF00003\t400\t500\n")
    return rows


@pytest.fixture(scope="module")
def dataset(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A working directory laid out the way the runner expects to find inputs."""
    root = tmp_path_factory.mktemp("dcgo_e2e")
    (root / "data" / "interim").mkdir(parents=True)
    (root / "data" / "raw" / "goa_annotations").mkdir(parents=True)

    gaf = root / "data" / "raw" / "goa_annotations" / f"goa_{SPECIES}.gaf.gz"
    with gzip.open(gaf, "wt") as handle:
        handle.writelines(_gaf_rows(N_POSITIVE, N_NEGATIVE))

    ipr = root / "data" / "interim" / f"protein2ipr_{SPECIES}.dat.gz"
    with gzip.open(ipr, "wt") as handle:
        handle.writelines(_protein2ipr_rows(N_POSITIVE, N_NEGATIVE))

    return root


def _cli_command() -> list[str]:
    """Prefer the installed `dcgo` console script; fall back to the script path.

    The P0 item is specifically about the *installed* command, so when the
    project is installed (the normal `uv sync` case, and what CI builds a wheel
    for) this exercises the same entry point a user gets. The fallback keeps the
    test runnable from a bare checkout.
    """
    console_script = Path(sys.executable).parent / "dcgo"
    if console_script.exists():
        return [str(console_script)]
    return [sys.executable, str(REPO_ROOT / "run_dcgo_human.py")]


@pytest.fixture(scope="module")
def completed_run(dataset: Path) -> tuple[subprocess.CompletedProcess, Path]:
    """Run the CLI once for the whole module; the assertions all read its output."""
    output_dir = dataset / "results"
    result = subprocess.run(
        _cli_command()
        + [
            "--species",
            SPECIES,
            "--output-dir",
            str(output_dir),
            "--num-cores",
            "1",
        ],
        cwd=dataset,
        capture_output=True,
        text=True,
        timeout=300,
        env={"PYTHONPATH": str(REPO_ROOT), "PATH": "/usr/bin:/bin"},
    )
    return result, output_dir


class TestCliRun:
    def test_exits_successfully(
        self, completed_run: tuple[subprocess.CompletedProcess, Path]
    ) -> None:
        result, _output_dir = completed_run
        assert result.returncode == 0, (
            f"CLI failed (exit {result.returncode})\n"
            f"--- stdout ---\n{result.stdout[-4000:]}\n"
            f"--- stderr ---\n{result.stderr[-4000:]}"
        )

    def test_writes_the_expected_files(
        self, completed_run: tuple[subprocess.CompletedProcess, Path]
    ) -> None:
        _result, output_dir = completed_run
        assert (output_dir / "domain_go_associations_significant.tsv").exists()
        assert (output_dir / "domain_go_associations_top100.tsv").exists()


class TestOutputSchema:
    def _rows(self, output_dir: Path) -> list[dict[str, str]]:
        path = output_dir / "domain_go_associations_significant.tsv"
        with open(path, newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_header_matches_the_documented_columns(
        self, completed_run: tuple[subprocess.CompletedProcess, Path]
    ) -> None:
        _result, output_dir = completed_run
        path = output_dir / "domain_go_associations_significant.tsv"
        header = path.read_text().splitlines()[0].split("\t")
        assert header == EXPECTED_COLUMNS

    def test_every_field_parses_as_its_type(
        self, completed_run: tuple[subprocess.CompletedProcess, Path]
    ) -> None:
        _result, output_dir = completed_run
        rows = self._rows(output_dir)
        assert rows, "no significant associations in a fixture built to have one"
        for row in rows:
            assert row["domain"].startswith("IPR")
            assert row["go_term"].startswith("GO:")
            assert 0.0 <= float(row["p_value"]) <= 1.0
            assert 0.0 <= float(row["adj_p_value"]) <= 1.0
            float(row["odds_ratio"])
            assert 0.0 <= float(row["hyper_score"]) <= 100.0
            assert row["domain_type"] in {"single", "supra_pair", "supra_triple"}
            assert int(row["n_observations"]) > 0

    def test_q_values_are_never_smaller_than_p_values(
        self, completed_run: tuple[subprocess.CompletedProcess, Path]
    ) -> None:
        """BH adjustment can only move a p-value up."""
        _result, output_dir = completed_run
        for row in self._rows(output_dir):
            assert float(row["adj_p_value"]) >= float(row["p_value"]) - 1e-12


class TestResultIsCorrect:
    """Schema-valid output is not enough; the answer has to be right."""

    def _pairs(self, output_dir: Path) -> set[tuple[str, str]]:
        path = output_dir / "domain_go_associations_significant.tsv"
        with open(path, newline="") as handle:
            return {
                (row["domain"], row["go_term"])
                for row in csv.DictReader(handle, delimiter="\t")
            }

    def test_finds_the_planted_association(
        self, completed_run: tuple[subprocess.CompletedProcess, Path]
    ) -> None:
        _result, output_dir = completed_run
        assert (IPR_HIT, TRUE_TERM) in self._pairs(output_dir)

    def test_does_not_report_the_uninformative_domain(
        self, completed_run: tuple[subprocess.CompletedProcess, Path]
    ) -> None:
        _result, output_dir = completed_run
        noise = {pair for pair in self._pairs(output_dir) if pair[0] == IPR_NOISE}
        assert not noise, f"uninformative domain reported as significant: {noise}"

    def test_supra_domains_reach_the_output(
        self, completed_run: tuple[subprocess.CompletedProcess, Path]
    ) -> None:
        """Supra-domains are on by default, so the output must contain some.

        The fixture gives every protein `IPR_UBIQUITOUS` next to whatever else
        it carries, so contiguous combinations exist to be generated. Without
        this the whole supra-domain half of the pipeline could break and the
        other assertions would still pass.
        """
        _result, output_dir = completed_run
        path = output_dir / "domain_go_associations_significant.tsv"
        with open(path, newline="") as handle:
            types = {
                row["domain_type"] for row in csv.DictReader(handle, delimiter="\t")
            }
        assert "single" in types
        assert types & {"supra_pair", "supra_triple"}, (
            f"no supra-domain associations in the output (types seen: {types})"
        )

    def test_the_planted_association_ranks_first(
        self, completed_run: tuple[subprocess.CompletedProcess, Path]
    ) -> None:
        _result, output_dir = completed_run
        path = output_dir / "domain_go_associations_top100.tsv"
        with open(path, newline="") as handle:
            first = next(csv.DictReader(handle, delimiter="\t"))
        assert (first["domain"], first["go_term"]) == (IPR_HIT, TRUE_TERM)


class TestFailsLoudly:
    """A missing input must be an error, not an empty-but-successful run."""

    def test_missing_species_input_is_a_nonzero_exit(self, tmp_path: Path) -> None:
        result = subprocess.run(
            _cli_command()
            + [
                "--species",
                "does-not-exist",
                "--output-dir",
                str(tmp_path / "results"),
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            timeout=120,
            env={"PYTHONPATH": str(REPO_ROOT), "PATH": "/usr/bin:/bin"},
        )
        assert result.returncode != 0
