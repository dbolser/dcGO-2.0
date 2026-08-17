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
        species="human",
        min_support=0,
        min_ic=0.0,
        ontology="go",
        enable_relative_inference=False,
        propagate_annotations=False,
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


class TestRejectsInvalidArguments:
    @pytest.mark.parametrize(
        "overrides, expected",
        [
            ({"fdr_threshold": 5.0}, "--fdr-threshold"),
            ({"fdr_threshold": 0.0}, "--fdr-threshold"),
            ({"fdr_threshold": -0.1}, "--fdr-threshold"),
            ({"batch_size": 0}, "--batch-size"),
            ({"batch_size": -1}, "--batch-size"),
            ({"min_support": -1}, "--min-support"),
            ({"min_ic": -0.5}, "--min-ic"),
            ({"num_cores": 0}, "--num-cores"),
            ({"num_cores": -4}, "--num-cores"),
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


class TestOddsRatioInterval:
    """The interval published beside each association's contingency cells."""

    def test_brackets_the_point_estimate(self) -> None:
        from run_dcgo_human import odds_ratio_interval

        a, b, c, d = 30, 5, 5, 200
        low, high = odds_ratio_interval(a, b, c, d)
        point = (a * d) / (b * c)
        assert low < point < high

    def test_a_sparse_table_gives_a_wide_interval(self) -> None:
        """The whole point: FDR significance alone hides fragility."""
        from run_dcgo_human import odds_ratio_interval

        sparse_low, sparse_high = odds_ratio_interval(2, 1, 1, 100)
        dense_low, dense_high = odds_ratio_interval(200, 100, 100, 10000)
        # Same odds ratio, two orders of magnitude apart in support.
        assert (sparse_high / sparse_low) > (dense_high / dense_low)

    def test_zero_cell_is_haldane_corrected_not_infinite(self) -> None:
        from run_dcgo_human import odds_ratio_interval

        low, high = odds_ratio_interval(10, 0, 5, 100)
        assert math.isfinite(low) and math.isfinite(high)
        assert low > 0

    def test_negative_cells_are_nan(self) -> None:
        from run_dcgo_human import odds_ratio_interval

        low, high = odds_ratio_interval(-1, 5, 5, 100)
        assert math.isnan(low) and math.isnan(high)

    def test_cells_reconstruct_the_reported_odds_ratio(self) -> None:
        """Publishing a/b/c/d is only useful if they reproduce the row.

        `fisher_exact_vectorized_batch` computes the odds ratio as (a*d)/(b*c);
        a reader recomputing from the exported cells must land on the same
        number, or the columns are decoration.
        """
        import numpy as np

        from src.vectorized_fisher import fisher_exact_vectorized_batch

        tables = np.array([[[30, 5], [5, 200]], [[7, 2], [3, 88]]], dtype=np.int32)
        odds_ratios, _pvalues = fisher_exact_vectorized_batch(tables)
        for table, reported in zip(tables, odds_ratios):
            a, b = int(table[0, 0]), int(table[0, 1])
            c, d = int(table[1, 0]), int(table[1, 1])
            assert (a * d) / (b * c) == pytest.approx(reported)
