"""Unit tests for the AnnotationSource seam (multi-ontology support)."""

import gzip

import pytest

from src.annotation_source import (
    GO_SPEC,
    AnnotationSource,
    GAFAnnotationSource,
    OntologySpec,
    restrict_to_universe,
)


@pytest.fixture
def sample_gaf_file(tmp_path):
    """A tiny GAF 2.2 file with one experimental and one IEA annotation."""
    file_path = tmp_path / "goa_mouse.gaf.gz"
    data = [
        "!gaf-version: 2.2",
        "UniProtKB\tP12345\tPROT1\t\tGO:0008150\tPMID:1\tIDA\t\tP\tProtein 1\t\tprotein\ttaxon:10090\t20210101\tGOA",
        "UniProtKB\tP12345\tPROT1\t\tGO:0005515\tPMID:2\tIEA\t\tF\tProtein 1\t\tprotein\ttaxon:10090\t20210101\tGOA",
        "UniProtKB\tP23456\tPROT2\t\tGO:0006810\tPMID:3\tIMP\t\tP\tProtein 2\t\tprotein\ttaxon:10090\t20210101\tGOA",
    ]
    with gzip.open(file_path, "wt") as f:
        f.write("\n".join(data))
    return file_path


class TestOntologySpec:
    def test_go_default_spec(self):
        assert GO_SPEC.ontology_id == "GO"
        assert GO_SPEC.term_prefix == "GO:"
        assert GO_SPEC.obo_path is None

    def test_custom_spec_without_prefix(self):
        # EC numbers carry no prefix, so term_prefix stays None.
        ec = OntologySpec(ontology_id="EC", name="Enzyme Commission")
        assert ec.term_prefix is None


class TestGAFAnnotationSource:
    def test_parse_returns_protein_term_map(self, sample_gaf_file):
        source = GAFAnnotationSource(sample_gaf_file, evidence_filter="manual")
        result = source.parse()

        # Manual filter keeps IDA/IMP, drops IEA.
        assert result["P12345"] == {"GO:0008150"}
        assert result["P23456"] == {"GO:0006810"}
        assert "GO:0005515" not in result["P12345"]

    def test_default_spec_is_go(self, sample_gaf_file):
        source = GAFAnnotationSource(sample_gaf_file)
        assert source.spec is GO_SPEC

    def test_is_an_annotation_source(self, sample_gaf_file):
        assert isinstance(GAFAnnotationSource(sample_gaf_file), AnnotationSource)


class TestCustomAnnotationSource:
    """A non-GO ontology plugs in by subclassing and returning a UniProt-keyed map."""

    def test_subclass_supplies_terms(self):
        class FakeECSource(AnnotationSource):
            spec = OntologySpec(ontology_id="EC", name="Enzyme Commission")

            def parse(self):
                return {"P12345": {"1.1.1.1"}, "P23456": {"2.7.7.7"}}

        source = FakeECSource()
        result = source.parse()
        assert isinstance(source, AnnotationSource)
        assert result["P12345"] == {"1.1.1.1"}
        assert source.spec.ontology_id == "EC"

    def test_cannot_instantiate_without_parse(self):
        class Incomplete(AnnotationSource):
            spec = OntologySpec(ontology_id="X", name="X")

        with pytest.raises(TypeError):
            Incomplete()


class TestRestrictToUniverse:
    """The Fisher universe is the intersection, not the union, of the two maps."""

    def test_drops_proteins_without_domains(self) -> None:
        annotations = {"P1": {"T1"}, "P2": {"T2"}, "P3": {"T3"}}
        restricted = restrict_to_universe(annotations, {"P1", "P3"})
        assert restricted == {"P1": {"T1"}, "P3": {"T3"}}

    def test_ignores_universe_members_with_no_annotations(self) -> None:
        restricted = restrict_to_universe({"P1": {"T1"}}, {"P1", "P2"})
        assert restricted == {"P1": {"T1"}}

    def test_excludes_other_species_from_the_fisher_universe(self) -> None:
        """The regression this exists for.

        ``uniprot_sprot.dat`` covers every organism while ``protein2ipr`` is
        extracted per species, so an unrestricted UniProt-native source put
        non-human proteins into every contingency table as domain-negative
        observations. Here HUMAN1/HUMAN2 have domains; MOUSE1/YEAST1 do not.
        """
        annotations = {
            "HUMAN1": {"R-HSA-1"},
            "HUMAN2": {"R-HSA-1"},
            "MOUSE1": {"R-HSA-1"},
            "YEAST1": {"R-HSA-2"},
        }
        domain_proteins = {"HUMAN1", "HUMAN2"}

        restricted = restrict_to_universe(annotations, domain_proteins)

        assert set(restricted) == domain_proteins
        # The term's background is now 2/2 human proteins rather than 3/4 across
        # three organisms, which is what the Fisher table needs.
        carriers = [p for p, terms in restricted.items() if "R-HSA-1" in terms]
        assert len(carriers) == 2

    def test_empty_universe_yields_empty_map(self) -> None:
        assert restrict_to_universe({"P1": {"T1"}}, set()) == {}
