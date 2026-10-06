"""Fixture-driven tests for validation programs that produce reported metrics."""

from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PARENT = "GO:0000001"
CHILD = "GO:0000002"
ROOT = "GO:0000003"
UNKNOWN = "GO:9999999"


def _run(script: str, *args: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(REPO_ROOT / "validation" / script), *map(str, args)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=60,
    )


@pytest.fixture
def go_obo(tmp_path: Path) -> Path:
    path = tmp_path / "go.obo"
    path.write_text(
        f"""format-version: 1.2

[Term]
id: {PARENT}
name: parent
namespace: biological_process

[Term]
id: {CHILD}
name: child
namespace: biological_process
is_a: {PARENT} ! parent

[Term]
id: {ROOT}
name: another root
namespace: biological_process
""",
        encoding="utf-8",
    )
    return path


def _gaf_row(protein: str, term: str) -> str:
    return "\t".join(
        [
            "UniProtKB",
            protein,
            protein,
            "",
            term,
            "PMID:1",
            "IDA",
            "",
            "P",
            "",
            "",
            "protein",
            "taxon:9606",
            "20200101",
            "GOA",
        ]
    )


def test_apply_relative_inference_filters_nonspecific_pair_and_keeps_untestable(
    tmp_path: Path, go_obo: Path
) -> None:
    predictions = tmp_path / "predictions.tsv"
    pd.DataFrame(
        [
            ("IPR_STRONG", CHILD, 0.001),
            ("IPR_WEAK", CHILD, 0.001),
            ("IPR_STRONG", ROOT, 0.002),
            ("IPR_STRONG", UNKNOWN, 0.003),
        ],
        columns=["domain", "go_term", "p_value"],
    ).to_csv(predictions, sep="\t", index=False)

    gaf = tmp_path / "t0.gaf.gz"
    with gzip.open(gaf, "wt") as handle:
        handle.write("!gaf-version: 2.2\n")
        for index in range(10):
            term = CHILD if index < 5 else PARENT
            handle.write(_gaf_row(f"P{index}", term) + "\n")

    interpro = tmp_path / "protein2ipr.dat.gz"
    with gzip.open(interpro, "wt") as handle:
        for index in range(10):
            if index < 5:
                handle.write(f"P{index}\tIPR_STRONG\tStrong\tPF00001\t10\t100\n")
            if index in {0, 1, 5, 6, 7}:
                handle.write(f"P{index}\tIPR_WEAK\tWeak\tPF00002\t120\t220\n")

    output = tmp_path / "relative.tsv"
    result = _run(
        "apply_relative_inference.py",
        "--predictions",
        predictions,
        "--t0-gaf",
        gaf,
        "--interpro",
        interpro,
        "--go-ontology",
        go_obo,
        "--output",
        output,
        "--min-background",
        10,
        "--disable-supra-domains",
    )

    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(output, sep="\t")
    assert set(zip(frame["domain"], frame["go_term"])) == {
        ("IPR_STRONG", CHILD),
        ("IPR_STRONG", ROOT),
        ("IPR_STRONG", UNKNOWN),
    }
    strong = frame[frame["go_term"] == CHILD].iloc[0]
    assert strong["overall_p"] == pytest.approx(0.001)
    assert strong["relative_p"] == pytest.approx(1 / 252)
    assert strong["p_value"] == pytest.approx(1 / 252)
    assert frame.loc[frame["go_term"].isin([ROOT, UNKNOWN]), "relative_p"].isna().all()


def test_domain_centric_cli_writes_true_path_aware_metrics(
    tmp_path: Path, go_obo: Path
) -> None:
    predictions = tmp_path / "predictions.tsv"
    pd.DataFrame(
        [
            ("IPR000001", CHILD),
            ("IPR000001", UNKNOWN),
            ("IPR000002", CHILD),
        ],
        columns=["domain", "go_term"],
    ).to_csv(predictions, sep="\t", index=False)
    reference = tmp_path / "interpro2go"
    reference.write_text(
        f"InterPro:IPR000001 First > GO:child ; {CHILD}\n"
        f"InterPro:IPR000002 Second > GO:parent ; {PARENT}\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "metrics"

    result = _run(
        "domain_centric_eval.py",
        "--predictions",
        f"fixture={predictions}",
        "--reference",
        reference,
        "--go-ontology",
        go_obo,
        "--output-dir",
        output_dir,
    )

    assert result.returncode == 0, result.stderr
    row = pd.read_csv(output_dir / "domain_centric_metrics.tsv", sep="\t").iloc[0]
    assert row["config"] == "fixture"
    assert row["n_associations"] == 3
    assert row["shared_domains"] == 2
    assert row["predicted_pairs"] == 5
    assert row["recovered"] == 3
    assert row["precision_lb"] == pytest.approx(3 / 5)
    assert row["recall"] == pytest.approx(1.0)
    assert row["f1"] == pytest.approx(0.75)


def test_domain_centric_cli_rejects_missing_association_columns(
    tmp_path: Path, go_obo: Path
) -> None:
    predictions = tmp_path / "bad.tsv"
    predictions.write_text("domain\tterm\nIPR000001\tGO:0000001\n", encoding="utf-8")
    reference = tmp_path / "interpro2go"
    reference.write_text(
        f"InterPro:IPR000001 First > GO:parent ; {PARENT}\n", encoding="utf-8"
    )

    result = _run(
        "domain_centric_eval.py",
        "--predictions",
        f"bad={predictions}",
        "--reference",
        reference,
        "--go-ontology",
        go_obo,
        "--output-dir",
        tmp_path / "output",
    )

    assert result.returncode == 1
    assert "missing required columns: go_term" in result.stderr


def test_validate_results_direct_script_uses_its_own_checkout(tmp_path: Path) -> None:
    """The script must not resolve ``validation`` from another checkout/PYTHONPATH."""
    foreign_root = tmp_path / "foreign-checkout"
    foreign_package = foreign_root / "validation"
    foreign_package.mkdir(parents=True)
    (foreign_package / "__init__.py").write_text("", encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "validation" / "validate_results.py")],
        cwd=tmp_path,
        env={"PYTHONPATH": str(foreign_root)},
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert result.returncode == 1
    assert "Missing required inputs" in result.stderr
    assert "ModuleNotFoundError" not in result.stderr


@pytest.fixture
def t0_benchmark(tmp_path: Path) -> dict[str, Path]:
    """A two-protein 2021→2026 split: P0 annotated at t0, P1 only at t1 (the
    held-out cohort), both carrying IPR_A in the protein2ipr subset."""
    t0 = tmp_path / "t0.gaf.gz"
    with gzip.open(t0, "wt") as handle:
        handle.write("!gaf-version: 2.2\n" + _gaf_row("P0", CHILD) + "\n")
    t1 = tmp_path / "t1.gaf.gz"
    with gzip.open(t1, "wt") as handle:
        handle.write("!gaf-version: 2.2\n")
        for protein in ("P0", "P1"):
            handle.write(_gaf_row(protein, CHILD) + "\n")
    interpro = tmp_path / "protein2ipr_t0.dat.gz"
    with gzip.open(interpro, "wt") as handle:
        for protein in ("P0", "P1"):
            handle.write(f"{protein}\tIPR_A\tA\tPF00001\t10\t100\n")
    return {"t0": t0, "t1": t1, "interpro": interpro}


def _write_run(directory: Path, interpro_sha256: str) -> None:
    """A rung's association tables and the manifest naming its protein2ipr."""
    directory.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(
        [("IPR_A", CHILD, 0.001, 0.01, 0.01)],
        columns=["domain", "go_term", "p_value", "adj_p_value", "q_value"],
    )
    for name in ("associations_significant", "annotations_propagated"):
        table.to_csv(directory / f"domain_go_{name}.tsv", sep="\t", index=False)
    (directory / "run_manifest_go.json").write_text(
        json.dumps(
            {
                "parameters": {"domain_key": "interpro", "evidence_filter": "manual"},
                "inputs": [{"role": "domain_annotations", "sha256": interpro_sha256}],
            }
        ),
        encoding="utf-8",
    )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_temporal_benchmark_scores_through_the_architectures_training_used(
    tmp_path: Path, go_obo: Path, t0_benchmark: dict[str, Path]
) -> None:
    run = tmp_path / "run"
    _write_run(run, _sha256(t0_benchmark["interpro"]))

    result = _run(
        "temporal_benchmark.py",
        "--t0-gaf",
        t0_benchmark["t0"],
        "--t1-gaf",
        t0_benchmark["t1"],
        "--predictions",
        run / "domain_go_associations_significant.tsv",
        "--interpro",
        t0_benchmark["interpro"],
        "--go-ontology",
        go_obo,
        "--n-permutations",
        1,
        "--min-ic",
        0,
        "--output-dir",
        tmp_path / "metrics",
    )

    assert result.returncode == 0, result.stderr
    metrics = pd.read_csv(
        tmp_path / "metrics" / "temporal_benchmark_metrics.tsv", sep="\t"
    )
    dcgo = metrics[metrics["method"] == "dcGO"].iloc[0]
    assert dcgo["n_eval_proteins"] == 1  # P1, transferred through IPR_A
    assert dcgo["f_max"] == pytest.approx(1.0)


def test_temporal_benchmark_refuses_other_architectures(
    tmp_path: Path, go_obo: Path, t0_benchmark: dict[str, Path]
) -> None:
    run = tmp_path / "run"
    _write_run(run, "0" * 64)  # trained on some other protein2ipr

    result = _run(
        "temporal_benchmark.py",
        "--t0-gaf",
        t0_benchmark["t0"],
        "--t1-gaf",
        t0_benchmark["t1"],
        "--predictions",
        run / "domain_go_associations_significant.tsv",
        "--interpro",
        t0_benchmark["interpro"],
        "--go-ontology",
        go_obo,
        "--output-dir",
        tmp_path / "metrics",
    )

    assert result.returncode == 1
    assert "protein2ipr sha256" in result.stderr
    assert not (tmp_path / "metrics").exists()


@pytest.mark.parametrize("matching", [True, False])
def test_ablation_scores_only_runs_trained_on_its_protein2ipr(
    tmp_path: Path, go_obo: Path, t0_benchmark: dict[str, Path], matching: bool
) -> None:
    sha256 = _sha256(t0_benchmark["interpro"])
    for rung in (
        "single",
        "supra",
        "supra_input",
        "supra_relative",
        "supra_output",
        "supra_input_relative",
        "supra_input_output",
        "supra_relative_output",
        "full",
    ):
        recorded = sha256 if matching or rung != "full" else "0" * 64
        _write_run(tmp_path / "runs" / rung, recorded)

    result = _run(
        "ablation.py",
        "--t0-gaf",
        t0_benchmark["t0"],
        "--t1-gaf",
        t0_benchmark["t1"],
        "--run-dir",
        tmp_path / "runs",
        "--interpro",
        t0_benchmark["interpro"],
        "--go-ontology",
        go_obo,
        "--min-ic",
        0,
        "--n-bootstrap",
        2,
        "--n-permutations",
        0,
        "--output-dir",
        tmp_path / "eval",
    )

    if not matching:
        assert result.returncode == 1
        assert "full: run protein2ipr sha256=" in result.stderr
        assert not (tmp_path / "eval" / "ablation_provenance.tsv").exists()
        return
    assert result.returncode == 0, result.stderr
    provenance = pd.read_csv(tmp_path / "eval" / "ablation_provenance.tsv", sep="\t")
    assert set(provenance["interpro_file"]) == {str(t0_benchmark["interpro"])}
    assert set(provenance["interpro_sha256"]) == {sha256}
