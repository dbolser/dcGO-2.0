"""Release-stamped filenames shared by config URLs and CLI flag defaults.

WormBase and FlyBase put the release in the *filename*
(``phenotype_association.WS298.wb.gz``, ``genotype_phenotype_data_fb_2026_02
.tsv.gz``), so the same string must appear in two places: the download URL in
``config/settings.py`` and the on-disk default of the matching
``run_dcgo_human.py`` flag. Deriving both from one constant here makes a
release bump a single edit per source — change the release string, and the
URL and the flag default move together.

This lives in ``src`` (not ``config``) because ``run_dcgo_human.py`` must not
import ``config.settings`` at module scope: config is a source-checkout
convenience with directory-creating side effects, while ``src`` is the shipped
package.
"""

from __future__ import annotations

#: WormBase release; downloads.wormbase.org serves the versioned filename
#: under current-production-release/ (the per-release archive paths 403).
WORMBASE_RELEASE = "WS298"
WORMBASE_PHENOTYPE_FILENAME = f"phenotype_association.{WORMBASE_RELEASE}.wb.gz"
WORMBASE_ANATOMY_FILENAME = f"anatomy_association.{WORMBASE_RELEASE}.wb.gz"

#: FlyBase release, as it appears in URL paths (``FB2026_02``); filenames
#: carry it as ``fb_2026_02`` (underscore after the prefix).
FLYBASE_RELEASE = "FB2026_02"
_FB = f"fb_{FLYBASE_RELEASE.removeprefix('FB').lower()}"
FLYBASE_GENOTYPE_PHENOTYPE_FILENAME = f"genotype_phenotype_data_{_FB}.tsv.gz"
FLYBASE_FBAL_TO_FBGN_FILENAME = f"fbal_to_fbgn_{_FB}.tsv.gz"
FLYBASE_FBGN_UNIPROT_FILENAME = f"fbgn_NAseq_Uniprot_{_FB}.tsv.gz"
