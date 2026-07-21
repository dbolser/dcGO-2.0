"""Unit tests for configuration helpers (per-species GOA URLs)."""

import pytest

from config.settings import ConfigurationError, config


class TestGoaUrlFor:
    def test_human_matches_pinned_source(self):
        # For human, the templated URL must reproduce the pinned data source.
        assert config.goa_url_for("human") == config.data_sources["goa_annotations"].url

    @pytest.mark.parametrize(
        "species,expected",
        [
            ("mouse", "MOUSE/goa_mouse.gaf.gz"),
            ("zebrafish", "ZEBRAFISH/goa_zebrafish.gaf.gz"),
            ("rat", "RAT/goa_rat.gaf.gz"),
        ],
    )
    def test_species_layout(self, species, expected):
        url = config.goa_url_for(species)
        assert url == f"{config.goa_base_url}/{expected}"

    def test_case_insensitive(self):
        assert config.goa_url_for("Mouse") == config.goa_url_for("mouse")

    def test_empty_species_rejected(self):
        with pytest.raises(ConfigurationError):
            config.goa_url_for("  ")
