"""The two hierarchy post-processing stages are selected independently.

`--enable-true-path` and `--enable-relative-inference` are two different steps
of the dcGO paper — Step 3 (propagate to ancestors, only ever *adds*
annotations) and Step 2's relative inference (a parental-background Fisher
test, only ever *removes* associations). They shared one flag until this split,
which is why the §4 ablation could not attribute its result to either.

These tests pin the CLI contract and the request translation. The numerical
behaviour of each stage is covered by `test_ontology_processor.py` (filter) and
`test_hierarchy.py` (propagation).
"""

from __future__ import annotations

import argparse

import pytest

from run_dcgo_human import build_argument_parser, parse_arguments, validate_arguments
from src.runner import parse_run_request


def _args(**overrides) -> argparse.Namespace:
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


class TestFlagsAreIndependent:
    def test_both_flags_exist_and_default_off(self) -> None:
        args, _ = parse_arguments([])
        assert args.enable_true_path is False
        assert args.enable_relative_inference is False

    @pytest.mark.parametrize(
        "argv, true_path, relative",
        [
            (["--enable-true-path"], True, False),
            (["--enable-relative-inference"], False, True),
            (["--enable-true-path", "--enable-relative-inference"], True, True),
        ],
    )
    def test_each_combination_is_reachable(
        self, argv: list[str], true_path: bool, relative: bool
    ) -> None:
        args, _ = parse_arguments(argv)
        assert args.enable_true_path is true_path
        assert args.enable_relative_inference is relative

    def test_true_path_alone_no_longer_implies_the_filter(self) -> None:
        """The behaviour change this split exists for.

        Before the split, `--enable-true-path` on `--ontology go` silently also
        ran the parental-background filter, discarding roughly half the
        associations. It now propagates and nothing else.
        """
        args, _ = parse_arguments(["--enable-true-path"])
        assert args.enable_relative_inference is False

    def test_help_states_that_true_path_is_propagation_only(self) -> None:
        help_text = build_argument_parser().format_help()
        assert "--enable-relative-inference" in help_text
        assert "only" in help_text


class TestRelativeInferenceNeedsAHierarchy:
    """The filter ranges over a term's *direct parents*, so it needs a hierarchy.

    It is no longer GO-only: every registry ontology supplying ``build_parents``
    can run it. Only the flat cross-reference layers are rejected.
    """

    @pytest.mark.parametrize(
        "ontology",
        ["go", "ec", "reactome", "keyword", "doid", "tcdb", "subcellular", "ligand"],
    )
    def test_hierarchical_ontologies_are_accepted(self, ontology: str) -> None:
        parser = argparse.ArgumentParser(prog="dcgo")
        validate_arguments(
            _args(ontology=ontology, enable_relative_inference=True), parser
        )

    @pytest.mark.parametrize("ontology", ["disease", "complex", "rhea", "xref"])
    def test_flat_vocabularies_are_rejected(self, ontology: str) -> None:
        parser = argparse.ArgumentParser(prog="dcgo")
        with pytest.raises(SystemExit) as excinfo:
            validate_arguments(
                _args(ontology=ontology, enable_relative_inference=True), parser
            )
        assert excinfo.value.code == 2

    def test_a_flat_vocabulary_without_the_flag_is_unaffected(self) -> None:
        parser = argparse.ArgumentParser(prog="dcgo")
        validate_arguments(
            _args(ontology="disease", enable_relative_inference=False), parser
        )

    def test_direct_parents_accompany_every_hierarchy(self) -> None:
        """``build_ancestors`` is not a substitute; pin the two to each other."""
        from src.ontology_registry import ONTOLOGIES

        for key, entry in ONTOLOGIES.items():
            assert entry.supports_true_path == entry.supports_relative_inference, key


class TestRunRequestCarriesBothStages:
    def test_request_records_the_two_flags_separately(self) -> None:
        request = parse_run_request(["--enable-relative-inference"])
        assert request.enable_relative_inference is True
        assert request.enable_true_path is False

    def test_relative_inference_alone_still_requires_the_hierarchy_input(self) -> None:
        """The filter reads the same OBO propagation does, so it is mandatory."""
        from src.runner import resolve_inputs

        request = parse_run_request(
            ["--enable-relative-inference", "--go-ontology", "does-not-exist.obo"]
        )
        assert "go_obo" in " ".join(resolve_inputs(request).missing_inputs)


class TestInputPropagationIsAThirdIndependentStage:
    """`--propagate-annotations` is the True Path Rule on the *input* map.

    An annotation to a child term implies its parents by definition, so this is
    not a method choice about the domain associations — it is the correct reading
    of the protein annotation data, applied before any test is built. It is
    distinct from `--enable-true-path`, which propagates the *inferred*
    associations afterwards.
    """

    def test_it_defaults_off_and_is_independent(self) -> None:
        args, _ = parse_arguments([])
        assert args.propagate_annotations is False

        args, _ = parse_arguments(["--propagate-annotations"])
        assert args.propagate_annotations is True
        assert args.enable_relative_inference is False
        assert args.enable_true_path is False

    @pytest.mark.parametrize(
        "ontology", ["go", "ec", "reactome", "doid", "subcellular", "ligand"]
    )
    def test_hierarchical_ontologies_are_accepted(self, ontology: str) -> None:
        parser = argparse.ArgumentParser(prog="dcgo")
        validate_arguments(_args(ontology=ontology, propagate_annotations=True), parser)

    @pytest.mark.parametrize("ontology", ["disease", "rhea", "xref"])
    def test_a_flat_vocabulary_is_rejected(self, ontology: str) -> None:
        """Without a hierarchy a child annotation implies nothing."""
        parser = argparse.ArgumentParser(prog="dcgo")
        with pytest.raises(SystemExit) as excinfo:
            validate_arguments(
                _args(ontology=ontology, propagate_annotations=True), parser
            )
        assert excinfo.value.code == 2

    def test_it_alone_requires_the_hierarchy_input(self) -> None:
        from src.runner import resolve_inputs

        request = parse_run_request(
            ["--propagate-annotations", "--go-ontology", "does-not-exist.obo"]
        )
        assert "go_obo" in " ".join(resolve_inputs(request).missing_inputs)

    def test_the_request_carries_all_three_stages_separately(self) -> None:
        request = parse_run_request(
            [
                "--propagate-annotations",
                "--enable-relative-inference",
                "--enable-true-path",
            ]
        )
        assert request.propagate_annotations is True
        assert request.enable_relative_inference is True
        assert request.enable_true_path is True
