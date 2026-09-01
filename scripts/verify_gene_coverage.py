"""Verify gene → UniProt mapping coverage for the gene-keyed annotation layers.

Parses each gene-keyed adapter against its acquisition inputs under
``data/raw`` and prints ``RemapCoverage.value_coverage`` — the fraction of
distinct source gene ids that reached at least one UniProt accession — plus
the layer's protein count and its most-used unmapped ids. These are the
figures cited by ``paper/EVIDENCE_LEDGER.md`` block L; every adapter logs the
same counts on a pipeline run, so a manifest-carrying run re-emits them.

Usage: ``uv run python scripts/verify_gene_coverage.py [data/raw]``
"""

import sys
from pathlib import Path
from typing import Dict, Optional, Protocol, Set

from src.remap import RemapCoverage

RAW = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/raw")
SPROT = RAW / "uniprot_sprot_dat" / "uniprot_sprot.dat.gz"


class GeneKeyedSource(Protocol):
    """The slice of a gene-keyed adapter this script reports on."""

    coverage: Optional[RemapCoverage]

    def parse(self) -> Dict[str, Set[str]]: ...


def report(name: str, source: GeneKeyedSource) -> None:
    ann = source.parse()
    cov = source.coverage
    if cov is None:
        raise RuntimeError(f"{name}: adapter recorded no remap coverage")
    print(
        f"{name}: {cov.n_mapped_values:,}/{cov.n_source_values:,} gene ids mapped "
        f"({100 * cov.value_coverage:.1f}%); proteins in layer: {len(ann):,}; "
        f"unmapped examples: {cov.unmapped_values[:5]}"
    )


def main() -> None:
    from src.mgi_annotation_source import MGIAnnotationSource
    from src.wormbase_annotation_source import WormBasePhenotypeAnnotationSource
    from src.zfin_annotation_source import ZFINAnatomyAnnotationSource
    from src.flybase_annotation_source import (
        FBBT_SPEC,
        FBCV_SPEC,
        FlyBasePhenotypeAnnotationSource,
    )
    from src.hpo_annotation_source import HPOAnnotationSource
    from src.syngo_annotation_source import SynGOAnnotationSource

    report(
        "mp",
        MGIAnnotationSource(
            RAW / "mgi_reports" / "MGI_GenePheno.rpt",
            RAW / "mgi_reports" / "MRK_SwissProt_TrEMBL.rpt",
        ),
    )
    report(
        "wbphenotype",
        WormBasePhenotypeAnnotationSource(
            RAW / "wormbase" / "phenotype_association.WS298.wb.gz",
            RAW / "uniprot_idmapping" / "CAEEL_6239_idmapping.dat.gz",
        ),
    )
    report(
        "zfa",
        ZFINAnatomyAnnotationSource(
            RAW / "zfin" / "phenoGeneCleanData_fish.txt",
            RAW / "zfin" / "uniprot.txt",
        ),
    )
    fb = (
        RAW / "flybase" / "genotype_phenotype_data_fb_2026_02.tsv.gz",
        RAW / "flybase" / "fbal_to_fbgn_fb_2026_02.tsv.gz",
        RAW / "flybase" / "fbgn_NAseq_Uniprot_fb_2026_02.tsv.gz",
    )
    report("fbcv", FlyBasePhenotypeAnnotationSource(*fb, spec=FBCV_SPEC))
    report("fbbt", FlyBasePhenotypeAnnotationSource(*fb, spec=FBBT_SPEC))
    report(
        "hpo",
        HPOAnnotationSource(RAW / "hpo" / "genes_to_phenotype.txt", SPROT),
    )
    report(
        "syngo",
        SynGOAnnotationSource(RAW / "syngo" / "syngo1.3_complete_data.zip", SPROT),
    )


if __name__ == "__main__":
    main()
