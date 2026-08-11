"""Public parser behavior and programmatic runner entry-point tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from run_dcgo_human import build_argument_parser, parse_arguments


def test_help_uses_the_real_parser_contract(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        parse_arguments(["--help"])

    assert excinfo.value.code == 0
    help_text = capsys.readouterr().out
    assert "--ontology NAME" in help_text
    assert "ontologies available to --ontology:" in help_text
    assert "--species SPECIES" in help_text


def test_invalid_choice_exits_two_and_names_the_option(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        parse_arguments(["--ontology", "not-an-ontology"])

    assert excinfo.value.code == 2
    assert "--ontology" in capsys.readouterr().err


def test_parse_arguments_can_be_called_without_process_global_argv() -> None:
    args, parser = parse_arguments(
        ["--ontology", "ec", "--fdr-threshold", "0.05", "--output-dir", "out"]
    )

    assert args.ontology == "ec"
    assert args.fdr_threshold == 0.05
    assert args.output_dir == Path("out")
    assert parser is not build_argument_parser()
