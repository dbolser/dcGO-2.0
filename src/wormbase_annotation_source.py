"""WormBase phenotype ontology: worm gene phenotypes, re-keyed to UniProt.

The dcGO trick applied to *C. elegans* (``--species worm --ontology
wbphenotype``): domain → WBPhenotype associations are learned on worm proteins,
and the species-agnostic domains carry them to any proteome.

Annotation input is WormBase's ``phenotype_association.<release>.wb.gz`` — a
GAF 2.0 file keyed by WBGene id (column 2) with the WBPhenotype term in the
GO-id column (5). Rows whose qualifier column carries ``NOT`` assert the
phenotype was *not* observed; they are negative evidence, not annotations, so
they are dropped and counted. The file is gene-level already (one row per gene
× phenotype × evidence), so no genotype policy is needed — WormBase has
collapsed the allele → gene step for us (the underlying variant is the
``With/From`` column, which this layer does not need).

The WBGene → UniProt translation is the per-organism idmapping file
(``CAEEL_6239_idmapping.dat.gz``, rows typed ``WormBase``), which covers
TrEMBL as well as Swiss-Prot — essential here, since most worm proteins are
unreviewed. Applied at parse time with the counted unmapped/one-to-many policy
(:func:`src.gene_mapping.remap_gene_annotations`).

The hierarchy is ``wbphenotype.obo``, read by the shared light OBO reader.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Set

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec
from src.gene_mapping import parse_idmapping_accession_map, remap_gene_annotations
from src.goa_parser import has_not_qualifier
from src.hierarchy import open_text
from src.remap import RemapCoverage

#: WormBase phenotype ontology terms, e.g. ``WBPhenotype:0000061``.
WBPHENOTYPE_SPEC = OntologySpec(
    ontology_id="WBPhenotype",
    name="WormBase Phenotype Ontology",
    term_prefix="WBPhenotype:",
)

# GAF column indices (0-based): DB object id, qualifier, term id.
_GENE_COL = 1
_QUALIFIER_COL = 3
_TERM_COL = 4


def parse_wb_phenotype_association(path: Path) -> Dict[str, Set[str]]:
    """Parse a WormBase ``phenotype_association`` GAF into ``{WBGene: {term}}``.

    ``NOT``-qualified rows are dropped and counted (negative evidence); any row
    whose term column is not a ``WBPhenotype:`` id is counted as malformed.
    """
    logger.info(f"Parsing WormBase phenotype annotations from {path}")
    gene_terms: Dict[str, Set[str]] = defaultdict(set)
    n_rows = 0
    n_not = 0
    n_malformed = 0
    with open_text(path, label="WormBase phenotype association file") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= _TERM_COL:
                n_malformed += 1
                continue
            gene = fields[_GENE_COL].strip()
            term = fields[_TERM_COL].strip()
            if not gene.startswith("WBGene") or not term.startswith("WBPhenotype:"):
                n_malformed += 1
                continue
            if has_not_qualifier(fields[_QUALIFIER_COL]):
                n_not += 1
                continue
            n_rows += 1
            gene_terms[gene].add(term)
    logger.info(
        f"  Rows kept: {n_rows:,}; NOT-qualified rows dropped: {n_not:,}; "
        f"malformed: {n_malformed:,}"
    )
    logger.info(
        f"  Genes: {len(gene_terms):,}; terms: "
        f"{len({term for terms in gene_terms.values() for term in terms}):,}"
    )
    return dict(gene_terms)


class WormBasePhenotypeAnnotationSource(AnnotationSource):
    """Domain annotations keyed by WormBase phenotype term.

    Reads the phenotype-association GAF and re-keys its WBGene ids to UniProt
    accessions (idmapping ``WormBase`` rows) before the statistics see them.
    """

    spec = WBPHENOTYPE_SPEC

    def __init__(self, association_path: Path, idmapping_path: Path) -> None:
        self.association_path = Path(association_path)
        self.idmapping_path = Path(idmapping_path)
        #: Populated by :meth:`parse`; its *values* range over WBGene ids, its
        #: *keys* over WBPhenotype terms (axis-swapped, see remap_gene_annotations).
        self.coverage: Optional[RemapCoverage] = None

    def parse(self) -> Dict[str, Set[str]]:
        gene_terms = parse_wb_phenotype_association(self.association_path)
        gene_map = parse_idmapping_accession_map(
            self.idmapping_path, "WormBase", id_space="WBGene"
        )
        remapped, self.coverage = remap_gene_annotations(
            gene_terms,
            gene_map,
            label="WBGene→UniProt (WBPhenotype)",
            target_label="UniProt accession",
        )
        return remapped
