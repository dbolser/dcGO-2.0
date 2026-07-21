"""Unit tests for the Enzyme Commission (EC) annotation source."""

import gzip
from types import SimpleNamespace

import pytest

from src.annotation_source import AnnotationSource
from src.ec_annotation_source import (
    EC_SPEC,
    ECAnnotationSource,
    ec_ancestors,
    parse_enzyme_dat,
    propagate_ec_annotations,
)


def _assoc(domain, ec, q_value, hyper_score=50.0):
    """A minimal stand-in for run_dcgo_human.AssociationResult."""
    return SimpleNamespace(
        domain=domain, go_term=ec, q_value=q_value, hyper_score=hyper_score
    )


# A miniature enzyme.dat covering: a normal multi-protein entry, a multi-line
# DR entry, a transferred entry (no DR), and a deleted entry (no DR).
SAMPLE_ENZYME_DAT = """CC   -----------------------------------------------------------------------
CC   Copyright notice header.
CC   -----------------------------------------------------------------------
//
ID   1.1.1.1
DE   Alcohol dehydrogenase.
AN   Aldehyde reductase.
DR   P07327, ADH1A_HUMAN;  P28469, ADH1A_MACMU;  Q5RBP7, ADH1A_PONAB;
DR   P25405, ADH1_ALLMI ;
//
ID   1.1.1.2
DE   Alcohol dehydrogenase (NADP(+)).
DR   P14550, AK1A1_HUMAN;
//
ID   1.1.1.3
DE   Transferred entry: 1.1.1.869, 1.1.1.870 and 1.1.1.871.
//
ID   1.1.1.74
DE   Deleted entry.
//
"""


@pytest.fixture
def enzyme_dat_file(tmp_path):
    path = tmp_path / "enzyme.dat"
    path.write_text(SAMPLE_ENZYME_DAT)
    return path


@pytest.fixture
def enzyme_dat_gz(tmp_path):
    path = tmp_path / "enzyme.dat.gz"
    with gzip.open(path, "wt") as f:
        f.write(SAMPLE_ENZYME_DAT)
    return path


class TestEcAncestors:
    def test_full_ec_number(self):
        assert ec_ancestors("1.1.1.1") == ["1.1.1.-", "1.1.-.-", "1.-.-.-"]

    def test_partial_ec_number(self):
        assert ec_ancestors("1.1.1.-") == ["1.1.-.-", "1.-.-.-"]

    def test_top_level_has_no_ancestors(self):
        assert ec_ancestors("1.-.-.-") == []

    def test_malformed_number(self):
        assert ec_ancestors("not.an.ec") == []
        assert ec_ancestors("1.1.1") == []


class TestParseEnzymeDat:
    def test_maps_accessions_to_ec(self, enzyme_dat_file):
        result = parse_enzyme_dat(enzyme_dat_file)
        assert result["P07327"] == {"1.1.1.1"}
        assert result["P14550"] == {"1.1.1.2"}

    def test_multiline_dr_all_captured(self, enzyme_dat_file):
        # All four accessions of the 1.1.1.1 entry, across two DR lines.
        result = parse_enzyme_dat(enzyme_dat_file)
        for acc in ("P07327", "P28469", "Q5RBP7", "P25405"):
            assert result[acc] == {"1.1.1.1"}

    def test_transferred_and_deleted_ignored(self, enzyme_dat_file):
        result = parse_enzyme_dat(enzyme_dat_file)
        all_ec = set().union(*result.values())
        assert "1.1.1.3" not in all_ec  # transferred
        assert "1.1.1.74" not in all_ec  # deleted

    def test_gzip_supported(self, enzyme_dat_gz):
        result = parse_enzyme_dat(enzyme_dat_gz)
        assert result["P07327"] == {"1.1.1.1"}

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_enzyme_dat(tmp_path / "nope.dat")


class TestECAnnotationSource:
    def test_is_annotation_source(self, enzyme_dat_file):
        source = ECAnnotationSource(enzyme_dat_file)
        assert isinstance(source, AnnotationSource)
        assert source.spec is EC_SPEC
        assert source.spec.ontology_id == "EC"

    def test_parse_returns_protein_ec_map(self, enzyme_dat_file):
        result = ECAnnotationSource(enzyme_dat_file).parse()
        assert result["P28469"] == {"1.1.1.1"}
        assert "P14550" in result


class TestPropagateECAnnotations:
    def test_direct_plus_ancestors(self):
        anns = propagate_ec_annotations([_assoc("D1", "1.1.1.1", 0.001)])
        pairs = {(a.domain, a.go_term, a.annotation_type) for a in anns}
        assert ("D1", "1.1.1.1", "direct") in pairs
        assert ("D1", "1.1.1.-", "propagated") in pairs
        assert ("D1", "1.1.-.-", "propagated") in pairs
        assert ("D1", "1.-.-.-", "propagated") in pairs
        # Every annotation records the specific EC number it came from.
        assert all(a.direct_source_term == "1.1.1.1" for a in anns)

    def test_no_duplicate_pairs(self):
        anns = propagate_ec_annotations([_assoc("D1", "1.1.1.1", 0.001)])
        keys = [(a.domain, a.go_term) for a in anns]
        assert len(keys) == len(set(keys))

    def test_shared_ancestor_attributed_to_most_significant(self):
        # Both EC numbers share ancestors 1.1.-.- and 1.-.-.-; the lower-q source wins.
        anns = propagate_ec_annotations(
            [
                _assoc("D1", "1.1.1.1", 0.05),
                _assoc("D1", "1.1.2.1", 0.001),
            ]
        )
        shared = [a for a in anns if a.go_term == "1.1.-.-"]
        assert len(shared) == 1
        assert shared[0].direct_source_term == "1.1.2.1"
        assert shared[0].annotation_type == "propagated"

    def test_domains_are_independent(self):
        anns = propagate_ec_annotations(
            [_assoc("D1", "1.1.1.1", 0.001), _assoc("D2", "1.1.1.1", 0.001)]
        )
        assert any(a.domain == "D1" and a.go_term == "1.1.1.-" for a in anns)
        assert any(a.domain == "D2" and a.go_term == "1.1.1.-" for a in anns)
