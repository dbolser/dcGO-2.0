"""Mammalian Phenotype Ontology: MGI mouse phenotypes, re-keyed to UniProt.

This is the dcGO trick applied to mouse: the domain → MP association is learned
on *mouse* proteins (``--species mouse --ontology mp``), and because domains
are species-agnostic the resulting associations annotate any protein carrying
the domain — including human ones.

Annotation input is MGI's ``MGI_GenePheno.rpt``: one row per genotype × MP
term, eight unnamed tab-separated columns::

    allelic composition | allele symbol(s) | allele id(s) | genetic background
    | MP id | PubMed | MGI marker accession id(s) | MGI genotype accession id

The marker column (7) is the gene. **Single-gene policy**: a phenotype observed
on a multi-gene genotype (pipe-separated markers, e.g. an ``Atm``/``Rad50``
double mutant) cannot be attributed to either gene alone, so those rows are
dropped — counted and logged, never silently. The file is already
overwhelmingly single-gene (283,001 of 283,003 rows at acquisition,
2 multi-gene rows dropped), because MGI publishes multi-gene genotypes in a
separate report.

The MGI → UniProt translation is ``MRK_SwissProt_TrEMBL.rpt`` (marker id in
column 1, space-separated Swiss-Prot/TrEMBL accessions in column 7), applied at
parse time with the DOID/HPO layers' counted policy for unmapped and
one-to-many ids (:func:`src.gene_mapping.remap_gene_annotations`): an MGI
marker cross-referenced by several accessions credits all of them, a marker
with none is dropped and counted.

The hierarchy is ``mp.obo``, read by the shared light OBO reader
(:func:`src.hierarchy.parse_obo_child_parents`).
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Set

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec
from src.gene_mapping import GeneAccessionMap, remap_gene_annotations
from src.remap import RemapCoverage

#: Mammalian Phenotype Ontology terms, e.g. ``MP:0001516``.
MP_SPEC = OntologySpec(
    ontology_id="MP", name="Mammalian Phenotype Ontology", term_prefix="MP:"
)

#: Column indices of MGI_GenePheno.rpt (the file has no header row).
_MP_TERM_COL = 4
_MARKER_COL = 6


def parse_mgi_genepheno(path: Path) -> Dict[str, Set[str]]:
    """Parse ``MGI_GenePheno.rpt`` into ``{MGI marker id: {MP term}}``.

    Multi-gene genotypes (pipe-separated marker column) are dropped and
    counted — see the module docstring for the policy.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"MGI GenePheno file not found: {path}")

    logger.info(f"Parsing MGI genotype→phenotype annotations from {path}")
    gene_terms: Dict[str, Set[str]] = defaultdict(set)
    n_rows = 0
    n_multi_gene = 0
    n_malformed = 0
    with open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= _MARKER_COL:
                n_malformed += 1
                continue
            marker = fields[_MARKER_COL].strip()
            term = fields[_MP_TERM_COL].strip()
            if not marker or not term.startswith("MP:"):
                n_malformed += 1
                continue
            if "|" in marker:
                n_multi_gene += 1
                continue
            n_rows += 1
            gene_terms[marker].add(term)
    logger.info(
        f"  Rows kept: {n_rows:,} (single-gene); multi-gene genotype rows "
        f"dropped: {n_multi_gene:,}; malformed: {n_malformed:,}"
    )
    logger.info(
        f"  Genes: {len(gene_terms):,}; terms: "
        f"{len({term for terms in gene_terms.values() for term in terms}):,}"
    )
    return dict(gene_terms)


def parse_mrk_swissprot(path: Path) -> GeneAccessionMap:
    """Parse ``MRK_SwissProt_TrEMBL.rpt`` into an MGI → accessions map.

    Column 1 is the MGI marker accession, column 7 the space-separated UniProt
    accessions (Swiss-Prot and TrEMBL). Markers listing several accessions are
    kept one-to-many, matching the counted expansion policy downstream.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"MGI SwissProt/TrEMBL report not found: {path}")

    logger.info(f"Building MGI → accession map from {path}")
    mapping: Dict[str, Set[str]] = defaultdict(set)
    n_rows = 0
    for line in open(path, "rt"):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7 or not fields[0].startswith("MGI:"):
            continue
        accessions = {acc for acc in fields[6].split() if acc}
        if accessions:
            n_rows += 1
            mapping[fields[0]] |= accessions
    gene_map = GeneAccessionMap("MGI", dict(mapping), n_rows)
    logger.info(
        f"  Markers with accessions: {len(gene_map):,} "
        f"({gene_map.n_one_to_many:,} one-to-many)"
    )
    return gene_map


class MGIAnnotationSource(AnnotationSource):
    """Domain annotations keyed by Mammalian Phenotype term.

    Reads ``MGI_GenePheno.rpt`` (single-gene genotypes only) and re-keys its
    MGI marker ids to UniProt accessions before the statistics see them.
    """

    spec = MP_SPEC

    def __init__(self, genepheno_path: Path, marker_map_path: Path) -> None:
        self.genepheno_path = Path(genepheno_path)
        self.marker_map_path = Path(marker_map_path)
        #: Populated by :meth:`parse`; its *values* range over MGI marker ids,
        #: its *keys* over MP terms (axis-swapped, see remap_gene_annotations).
        self.coverage: Optional[RemapCoverage] = None

    def parse(self) -> Dict[str, Set[str]]:
        gene_terms = parse_mgi_genepheno(self.genepheno_path)
        gene_map = parse_mrk_swissprot(self.marker_map_path)
        remapped, self.coverage = remap_gene_annotations(
            gene_terms, gene_map, label="MGI→UniProt (MP)"
        )
        return remapped
