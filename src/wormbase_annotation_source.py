"""WormBase gene → term layers, re-keyed to UniProt: phenotype and anatomy.

The dcGO trick applied to *C. elegans* (``--species worm``): domain → term
associations are learned on worm proteins, and the species-agnostic domains
carry them to any proteome. Two WormBase association files share one GAF 2.0
shape (WBGene id in column 2, term in the GO-id column 5) and one parser:

* ``phenotype_association.<release>.wb.gz`` (``--ontology wbphenotype``) —
  WBPhenotype terms, *phenotype* semantics ("mutating this gene perturbs
  this"). Rows whose qualifier column carries ``NOT`` assert the phenotype was
  not observed; negative evidence, dropped and counted.
* ``anatomy_association.<release>.wb.gz`` (``--ontology wbbt``) — WBbt anatomy
  terms, **expression** semantics: "expressed in this anatomical structure"
  (the file derives from WormBase's expression-pattern curation), like the HPA
  cell-type layer and unlike every phenotype layer. Its qualifier column
  grades the call — ``Certain`` / ``Enriched`` / ``Partial`` / ``Uncertain``
  / blank (52,744 / 376,854 / 3,411 / 4,217 / 23,997 rows at WS298).
  ``Uncertain`` rows are dropped and counted: an annotation the curators
  themselves flag as doubtful should not enter a Fisher table. The other
  grades all assert expression and are kept.

Both files are gene-level already (the underlying variant/pattern is the
``With/From`` column, which these layers do not need), so no genotype policy
applies — WormBase has collapsed the allele → gene step for us.

The WBGene → UniProt translation is the per-organism idmapping file
(``CAEEL_6239_idmapping.dat.gz``, rows typed ``WormBase``), which covers
TrEMBL as well as Swiss-Prot — essential here, since most worm proteins are
unreviewed. Applied at parse time with the counted unmapped/one-to-many policy
(:func:`src.gene_mapping.remap_gene_annotations`).

The hierarchies are ``wbphenotype.obo`` and ``wbbt.obo``, read by the shared
light OBO reader.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, FrozenSet, Optional, Set

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

#: WormBase anatomy ontology terms, e.g. ``WBbt:0003679``.
WBBT_SPEC = OntologySpec(
    ontology_id="WBbt",
    name="C. elegans Gross Anatomy Ontology",
    term_prefix="WBbt:",
)

# GAF column indices (0-based): DB object id, qualifier, term id.
_GENE_COL = 1
_QUALIFIER_COL = 3
_TERM_COL = 4


def _parse_wb_gaf(
    path: Path,
    term_prefix: str,
    label: str,
    dropped_qualifiers: FrozenSet[str] = frozenset(),
) -> Dict[str, Set[str]]:
    """Parse a WormBase GAF into ``{WBGene: {term}}``, with qualifier policy.

    ``NOT``-qualified rows are always dropped and counted (negative evidence);
    ``dropped_qualifiers`` names further qualifier values to drop (the anatomy
    file's ``Uncertain``). Rows whose term column does not carry
    ``term_prefix`` are counted as malformed.
    """
    logger.info(f"Parsing {label} from {path}")
    gene_terms: Dict[str, Set[str]] = defaultdict(set)
    n_rows = 0
    n_not = 0
    n_dropped_qualifier = 0
    n_malformed = 0
    with open_text(path, label=label) as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= _TERM_COL:
                n_malformed += 1
                continue
            gene = fields[_GENE_COL].strip()
            term = fields[_TERM_COL].strip()
            if not gene.startswith("WBGene") or not term.startswith(term_prefix):
                n_malformed += 1
                continue
            qualifier = fields[_QUALIFIER_COL].strip()
            if has_not_qualifier(qualifier):
                n_not += 1
                continue
            # GAF qualifiers are |-separated lists; match against the tokens,
            # not the whole field, so "NOT|Uncertain"-style compounds drop.
            if dropped_qualifiers and dropped_qualifiers & set(qualifier.split("|")):
                n_dropped_qualifier += 1
                continue
            n_rows += 1
            gene_terms[gene].add(term)
    logger.info(
        f"  Rows kept: {n_rows:,}; NOT-qualified rows dropped: {n_not:,}; "
        f"dropped qualifiers {sorted(dropped_qualifiers)}: "
        f"{n_dropped_qualifier:,}; malformed: {n_malformed:,}"
    )
    logger.info(
        f"  Genes: {len(gene_terms):,}; terms: "
        f"{len({term for terms in gene_terms.values() for term in terms}):,}"
    )
    return dict(gene_terms)


def parse_wb_phenotype_association(path: Path) -> Dict[str, Set[str]]:
    """Parse a WormBase ``phenotype_association`` GAF into ``{WBGene: {term}}``.

    ``NOT``-qualified rows are dropped and counted (negative evidence); any row
    whose term column is not a ``WBPhenotype:`` id is counted as malformed.
    """
    return _parse_wb_gaf(path, "WBPhenotype:", "WormBase phenotype association file")


def parse_wb_anatomy_association(path: Path) -> Dict[str, Set[str]]:
    """Parse a WormBase ``anatomy_association`` GAF into ``{WBGene: {term}}``.

    Expression semantics — see the module docstring. ``Uncertain``-qualified
    rows are dropped and counted alongside the (never yet observed in this
    file) ``NOT`` rows; ``Certain``/``Enriched``/``Partial``/blank are kept.
    """
    return _parse_wb_gaf(
        path,
        "WBbt:",
        "WormBase anatomy association file",
        dropped_qualifiers=frozenset({"Uncertain"}),
    )


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


class WormBaseAnatomyAnnotationSource(AnnotationSource):
    """Domain annotations keyed by WormBase anatomy term (expression-based).

    Reads the anatomy-association GAF — "expressed in", not "phenotype of";
    see the module docstring — and re-keys its WBGene ids to UniProt
    accessions (idmapping ``WormBase`` rows) before the statistics see them.
    """

    spec = WBBT_SPEC

    def __init__(self, association_path: Path, idmapping_path: Path) -> None:
        self.association_path = Path(association_path)
        self.idmapping_path = Path(idmapping_path)
        #: Populated by :meth:`parse`; its *values* range over WBGene ids, its
        #: *keys* over WBbt terms (axis-swapped, see remap_gene_annotations).
        self.coverage: Optional[RemapCoverage] = None

    def parse(self) -> Dict[str, Set[str]]:
        gene_terms = parse_wb_anatomy_association(self.association_path)
        gene_map = parse_idmapping_accession_map(
            self.idmapping_path, "WormBase", id_space="WBGene"
        )
        remapped, self.coverage = remap_gene_annotations(
            gene_terms,
            gene_map,
            label="WBGene→UniProt (WBbt)",
            target_label="UniProt accession",
        )
        return remapped
