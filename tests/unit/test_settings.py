"""Unit tests for configuration helpers (per-species GOA URLs, dataset checksums)."""

import pytest

from config.settings import (
    ALL_SPECIES_ALIASES,
    ConfigurationError,
    DataSource,
    config,
)


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

    @pytest.mark.parametrize("alias", sorted(ALL_SPECIES_ALIASES))
    def test_all_species_resolves_to_the_cross_organism_release(self, alias):
        # The multi-species background is not an organism: EBI publishes it as
        # one cross-organism file, not under an ALLSPECIES/ directory. The
        # per-species pattern would compose a plausible URL that 404s, and the
        # run manifest would then record it as the input's origin.
        assert (
            config.goa_url_for(alias)
            == f"{config.goa_base_url}/UNIPROT/goa_uniprot_all.gaf.gz"
        )

    def test_all_species_aliases_do_not_shadow_an_organism(self):
        # Guard against an alias later colliding with a real GOA species dir.
        assert not ALL_SPECIES_ALIASES & {"human", "mouse", "rat", "zebrafish"}


def _source(**kwargs) -> DataSource:
    return DataSource(name="t", url="https://example.org/f", description="d", **kwargs)


class TestChecksums:
    def test_none_when_unset(self):
        assert _source().checksum_parts() is None

    def test_algorithm_and_digest_split(self):
        assert _source(checksum="sha256:ABCD").checksum_parts() == ("sha256", "abcd")

    def test_bare_digest_defaults_to_sha256(self):
        assert _source(checksum="abcd").checksum_parts() == ("sha256", "abcd")

    def test_unknown_algorithm_rejected_at_construction(self):
        with pytest.raises(ConfigurationError, match="algorithm"):
            _source(checksum="crc32:abcd")

    def test_non_hex_digest_rejected_at_construction(self):
        # Better to fail here than after a multi-GB download.
        with pytest.raises(ConfigurationError, match="hex"):
            _source(checksum="sha256:not-a-digest")

    def test_disease_ontology_is_pinned_and_checksummed(self):
        # The reviewer's requirement: the DO release a run depends on must be
        # reproducible by content, not just by URL.
        source = config.data_sources["disease_ontology"]
        assert "/releases/" in source.url
        assert source.url.startswith("https://")
        algorithm, digest = source.checksum_parts()
        assert (algorithm, len(digest)) == ("sha256", 64)
