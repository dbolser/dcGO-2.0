"""Unit tests for the protein2ipr extract provenance sidecars."""

import json

import pytest

from src.universe_provenance import (
    ProvenanceConflictError,
    ensure_overwrite_allowed,
    marker_path,
    read_marker,
    write_marker,
)


@pytest.fixture
def extract(tmp_path):
    path = tmp_path / "protein2ipr_worm.dat.gz"
    path.write_bytes(b"fake gzip payload")
    return path


def write_goa_marker(extract):
    return write_marker(
        extract,
        selection_rule="goa",
        selection_sources=["data/raw/goa_annotations/goa_worm.gaf.gz"],
        interpro_source="data/raw/interpro_mappings/protein2ipr.dat.gz",
        n_accessions=100,
        n_matched_lines=4200,
        tool="extract_human_interpro.py",
        created="2026-08-17T00:00:00+00:00",
    )


class TestMarkerRoundTrip:
    def test_marker_sits_next_to_the_extract(self, extract):
        assert marker_path(extract).name == "protein2ipr_worm.dat.gz.provenance.json"

    def test_written_marker_reads_back(self, extract):
        write_goa_marker(extract)
        marker = read_marker(extract)
        assert marker.selection_rule == "goa"
        assert marker.n_matched_lines == 4200
        assert marker.created == "2026-08-17T00:00:00+00:00"

    def test_missing_marker_reads_as_none(self, extract):
        assert read_marker(extract) is None

    def test_unreadable_marker_reads_as_none(self, extract):
        marker_path(extract).write_text("not json")
        assert read_marker(extract) is None

    def test_alien_marker_reads_as_none(self, extract):
        marker_path(extract).write_text(json.dumps({"something": "else"}))
        assert read_marker(extract) is None


class TestOverwritePolicy:
    def test_fresh_path_is_allowed(self, tmp_path):
        ensure_overwrite_allowed(tmp_path / "new.dat.gz", "idmapping")

    def test_same_rule_refresh_is_allowed(self, extract):
        write_goa_marker(extract)
        ensure_overwrite_allowed(extract, "goa")

    def test_different_rule_is_refused(self, extract):
        write_goa_marker(extract)
        with pytest.raises(ProvenanceConflictError, match="--force"):
            ensure_overwrite_allowed(extract, "idmapping")

    def test_refusal_names_the_existing_rule_and_sources(self, extract):
        write_goa_marker(extract)
        with pytest.raises(ProvenanceConflictError, match="'goa'.*goa_worm.gaf.gz"):
            ensure_overwrite_allowed(extract, "idmapping")

    def test_force_overrides_the_refusal(self, extract):
        write_goa_marker(extract)
        ensure_overwrite_allowed(extract, "idmapping", force=True)

    def test_pre_marker_extract_is_allowed_with_warning(self, extract):
        # Legacy extracts predate the markers; refusing would strand them all.
        ensure_overwrite_allowed(extract, "idmapping")
