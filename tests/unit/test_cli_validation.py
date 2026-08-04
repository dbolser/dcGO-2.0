"""The runner rejects nonsensical parameters instead of running on them.

`ENGINEERING_SCIENTIFIC_REVIEW_TODOS.md` (P0): "Validate CLI parameters: FDR
range, positive batch size and core count, shrinkage range, supported
input/species selection …". argparse checks types, not ranges, so before this
an `--fdr-threshold 5` accepted every test as significant and a
`--shrinkage-strength 2` extrapolated past the prior — both of them silently,
and an hour into a run.

Also covers `calculate_hypergeometric_score`'s failure value. It used to return
50.0, which is a plausible medium-confidence score in the same 1–100 range the
function reports genuine results on, so a numerical failure was
indistinguishable from a real answer in the exported column.
"""

from __future__ import annotations

import argparse
import math

import pytest

from run_dcgo_human import calculate_hypergeometric_score, validate_arguments


def _args(**overrides) -> argparse.Namespace:
    """A parameter set that passes validation, with fields overridable."""
    defaults = dict(
        fdr_threshold=0.01,
        batch_size=50000,
        num_cores=8,
        shrinkage_strength=0.5,
        species="human",
    )
    defaults.update(overrides)
    return argparse.Namespace(**defaults)


@pytest.fixture
def parser() -> argparse.ArgumentParser:
    """A parser whose .error() raises SystemExit(2), as argparse's does."""
    return argparse.ArgumentParser(prog="dcgo")


class TestAcceptsValidArguments:
    def test_defaults_pass(self, parser: argparse.ArgumentParser) -> None:
        validate_arguments(_args(), parser)

    @pytest.mark.parametrize("threshold", [1e-12, 0.01, 0.05, 1.0])
    def test_fdr_boundaries(
        self, threshold: float, parser: argparse.ArgumentParser
    ) -> None:
        validate_arguments(_args(fdr_threshold=threshold), parser)

    @pytest.mark.parametrize("strength", [0.0, 0.5, 1.0])
    def test_shrinkage_boundaries(
        self, strength: float, parser: argparse.ArgumentParser
    ) -> None:
        """0 and 1 are both meaningful: no shrinkage, and full shrinkage."""
        validate_arguments(_args(shrinkage_strength=strength), parser)


class TestRejectsInvalidArguments:
    @pytest.mark.parametrize(
        "overrides, expected",
        [
            ({"fdr_threshold": 5.0}, "--fdr-threshold"),
            ({"fdr_threshold": 0.0}, "--fdr-threshold"),
            ({"fdr_threshold": -0.1}, "--fdr-threshold"),
            ({"batch_size": 0}, "--batch-size"),
            ({"batch_size": -1}, "--batch-size"),
            ({"num_cores": 0}, "--num-cores"),
            ({"num_cores": -4}, "--num-cores"),
            ({"shrinkage_strength": 1.5}, "--shrinkage-strength"),
            ({"shrinkage_strength": -0.5}, "--shrinkage-strength"),
            ({"species": ""}, "--species"),
            ({"species": "../etc/passwd"}, "--species"),
        ],
    )
    def test_exits_two_and_names_the_option(
        self,
        overrides: dict,
        expected: str,
        parser: argparse.ArgumentParser,
        capsys: pytest.CaptureFixture,
    ) -> None:
        with pytest.raises(SystemExit) as excinfo:
            validate_arguments(_args(**overrides), parser)
        assert excinfo.value.code == 2
        assert expected in capsys.readouterr().err

    def test_fdr_above_one_explains_why(
        self, parser: argparse.ArgumentParser, capsys: pytest.CaptureFixture
    ) -> None:
        """The message should say what goes wrong, not just that it is invalid."""
        with pytest.raises(SystemExit):
            validate_arguments(_args(fdr_threshold=5.0), parser)
        assert "every test significant" in capsys.readouterr().err


class TestHypergeometricScoreFailureValue:
    def test_returns_nan_not_a_plausible_score(self) -> None:
        """A failure must be visibly missing, not mid-range.

        Negative cell counts are not a reachable contingency table; they stand
        in for whatever numerical failure trips the except branch.
        """
        score = calculate_hypergeometric_score(-1, -1, -1, -1)
        assert math.isnan(score), f"expected NaN for an unusable table, got {score}"

    def test_empty_margins_are_zero_not_nan(self) -> None:
        """A legitimately empty overlap is a real answer of zero, not a failure."""
        assert calculate_hypergeometric_score(0, 10, 10, 100) == 0.0

    def test_a_real_table_scores_in_range(self) -> None:
        score = calculate_hypergeometric_score(30, 5, 5, 200)
        assert 1.0 <= score <= 100.0
