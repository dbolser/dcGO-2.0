"""The release pins must produce the exact upstream filenames.

These strings are load-bearing twice over — they are the download URLs' tails
and the CLI flag defaults — and the FlyBase naming is a trap (path segment
``FB2026_02``, filename infix ``fb_2026_02``), so the derivations are asserted
literally.
"""

from src.release_pins import (
    FLYBASE_FBAL_TO_FBGN_FILENAME,
    FLYBASE_FBGN_UNIPROT_FILENAME,
    FLYBASE_GENOTYPE_PHENOTYPE_FILENAME,
    FLYBASE_RELEASE,
    WORMBASE_PHENOTYPE_FILENAME,
    WORMBASE_RELEASE,
)


def test_wormbase_filename_carries_the_release():
    assert (
        WORMBASE_PHENOTYPE_FILENAME == f"phenotype_association.{WORMBASE_RELEASE}.wb.gz"
    )


def test_flybase_filenames_use_the_underscored_lowercase_infix():
    # FB2026_02 (path) vs fb_2026_02 (filename) — the mismatch that bit once.
    assert FLYBASE_RELEASE == "FB2026_02"
    assert FLYBASE_GENOTYPE_PHENOTYPE_FILENAME == (
        "genotype_phenotype_data_fb_2026_02.tsv.gz"
    )
    assert FLYBASE_FBAL_TO_FBGN_FILENAME == "fbal_to_fbgn_fb_2026_02.tsv.gz"
    assert FLYBASE_FBGN_UNIPROT_FILENAME == "fbgn_NAseq_Uniprot_fb_2026_02.tsv.gz"


def test_settings_urls_end_with_the_pinned_filenames():
    from config.settings import Config

    sources = Config(use_env_overrides=False).data_sources
    assert sources["wormbase_phenotype"].url.endswith(WORMBASE_PHENOTYPE_FILENAME)
    assert sources["flybase_genotype_phenotype"].url.endswith(
        FLYBASE_GENOTYPE_PHENOTYPE_FILENAME
    )
    assert sources["flybase_fbal_to_fbgn"].url.endswith(FLYBASE_FBAL_TO_FBGN_FILENAME)
    assert sources["flybase_fbgn_uniprot"].url.endswith(FLYBASE_FBGN_UNIPROT_FILENAME)
    # And the release-stamped URLs carry a staleness hint for the downloader.
    assert sources["wormbase_phenotype"].update_hint
