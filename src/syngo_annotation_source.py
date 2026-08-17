"""SynGO: expert-curated synaptic gene annotations, re-keyed to UniProt.

SynGO (Koopmans et al. 2019; syngoportal.org) curates synapse-specific
annotations against an ontology of GO cellular-component/biological-process
terms extended with SynGO-specific ``SYNGO:`` terms. The bulk release is one
zip of three xlsx sheets:

* ``annotations.xlsx`` — one row per curated annotation, keyed by HGNC id and
  symbol. The sheet's ``uniprot_id`` column is deliberately **not** used: it
  names the protein of the underlying *evidence experiment*, which is often a
  mouse or rat orthologue — keying on it would leak non-human accessions into
  a human universe. The HGNC gene is the curated subject, so it is what gets
  mapped (HGNC id first, approved-symbol fallback) to UniProt via the
  Swiss-Prot flat file (:mod:`src.gene_mapping`).
* ``ontologies.xlsx`` — every SynGO term with its ``parent_id``, i.e. the
  child→parent hierarchy ships in the same zip; no OBO needed.
* ``genes.xlsx`` — the curated gene universe (not needed here).

Terms keep their released ids (``GO:``/``SYNGO:`` mixed), so no term mapping
is involved — only the gene→accession re-key, with the standard counted policy
for unmapped and one-to-many ids.
"""

from __future__ import annotations

import io
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec
from src.gene_mapping import (
    GeneAccessionIndex,
    GeneAccessionMap,
    parse_gene_accession_index,
    remap_gene_annotations,
)
from src.remap import RemapCoverage

#: SynGO terms are GO ids plus SynGO-specific ``SYNGO:`` extensions, so the
#: spec declares no single term prefix.
SYNGO_SPEC = OntologySpec(ontology_id="SynGO", name="SynGO synaptic gene ontology")

_ANNOTATIONS_MEMBER = "annotations.xlsx"
_ONTOLOGIES_MEMBER = "ontologies.xlsx"


def _iter_sheet_rows(zip_path: Path, member: str) -> Iterator[Dict[str, object]]:
    """Yield each data row of one xlsx inside the SynGO zip, keyed by header.

    openpyxl needs a seekable stream, so the member is read into memory —
    the whole zip is under a megabyte.
    """
    import openpyxl

    zip_path = Path(zip_path)
    if not zip_path.exists():
        raise FileNotFoundError(f"SynGO release zip not found: {zip_path}")
    with zipfile.ZipFile(zip_path) as archive:
        payload = archive.read(member)
    workbook = openpyxl.load_workbook(io.BytesIO(payload), read_only=True)
    rows = workbook.active.iter_rows(values_only=True)
    header: List[str] = [str(cell) for cell in next(rows)]
    for row in rows:
        yield dict(zip(header, row))


def _cell(row: Dict[str, object], column: str) -> str:
    value = row.get(column)
    return str(value).strip() if value is not None else ""


def parse_syngo_annotations(
    zip_path: Path,
) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
    """Parse ``annotations.xlsx`` into ``({hgnc id: {term}}, {hgnc id: symbol})``.

    A few annotation rows name several genes in one cell
    (``HGNC:243;HGNC:244`` — paralog pairs the evidence could not separate);
    each listed gene gets the term. Symbols are paired positionally with the
    ids when the two cells split into the same number of parts, and feed the
    approved-symbol fallback for HGNC ids the flat file does not
    cross-reference.
    """
    logger.info(f"Parsing SynGO annotations from {zip_path}")
    gene_terms: Dict[str, Set[str]] = defaultdict(set)
    symbols: Dict[str, str] = {}
    n_rows = 0
    for row in _iter_sheet_rows(zip_path, _ANNOTATIONS_MEMBER):
        term = _cell(row, "go_id")
        gene_ids = [
            part.strip() for part in _cell(row, "hgnc_id").split(";") if part.strip()
        ]
        if not gene_ids or not term:
            continue
        n_rows += 1
        row_symbols = [part.strip() for part in _cell(row, "hgnc_symbol").split(";")]
        for position, gene_id in enumerate(gene_ids):
            gene_terms[gene_id].add(term)
            if len(row_symbols) == len(gene_ids) and row_symbols[position]:
                symbols[gene_id] = row_symbols[position]
    logger.info(
        f"  Rows: {n_rows:,}; genes: {len(gene_terms):,}; terms: "
        f"{len({term for terms in gene_terms.values() for term in terms}):,}"
    )
    return dict(gene_terms), symbols


def parse_syngo_hierarchy(zip_path: Path) -> Dict[str, Set[str]]:
    """Read ``ontologies.xlsx`` into a ``{child: {parents}}`` map.

    Each term row carries a single ``parent_id`` (the release is a tree with
    two roots — the synapse CC term and the synaptic-process BP term, whose
    rows have no parent).
    """
    child_to_parents: Dict[str, Set[str]] = {}
    for row in _iter_sheet_rows(zip_path, _ONTOLOGIES_MEMBER):
        term = _cell(row, "id")
        parent = _cell(row, "parent_id")
        if term and parent:
            child_to_parents.setdefault(term, set()).add(parent)
    return child_to_parents


def resolve_hgnc_accessions(
    gene_terms: Dict[str, Set[str]],
    symbols: Dict[str, str],
    index: GeneAccessionIndex,
) -> Tuple[GeneAccessionMap, int]:
    """Resolve each annotated HGNC id to accessions, symbol as fallback.

    Returns the composite :class:`GeneAccessionMap` restricted to the annotated
    genes, plus how many of them only resolved through their symbol.
    """
    resolved: Dict[str, Set[str]] = {}
    n_symbol_fallback = 0
    for gene in gene_terms:
        accessions = index.hgnc.targets(gene)
        if not accessions:
            symbol = symbols.get(gene, "")
            accessions = index.symbol.targets(symbol) if symbol else set()
            if accessions:
                n_symbol_fallback += 1
        if accessions:
            resolved[gene] = accessions
    if n_symbol_fallback:
        logger.info(
            f"  {n_symbol_fallback:,} HGNC ids resolved via their approved "
            "symbol (no DR HGNC id match)"
        )
    return (
        GeneAccessionMap("HGNC(+symbol)", resolved, index.hgnc.n_entries),
        n_symbol_fallback,
    )


class SynGOAnnotationSource(AnnotationSource):
    """Domain annotations keyed by SynGO term.

    Reads the SynGO bulk-release zip and re-keys its HGNC genes to UniProt
    accessions before the annotations reach the statistics.
    """

    spec = SYNGO_SPEC

    def __init__(self, zip_path: Path, dat_path: Path) -> None:
        self.zip_path = Path(zip_path)
        self.dat_path = Path(dat_path)
        #: Populated by :meth:`parse`; its *values* are the HGNC ids, its
        #: *keys* the SynGO terms (see :func:`src.gene_mapping.remap_gene_annotations`).
        self.coverage: Optional[RemapCoverage] = None
        #: HGNC ids that only resolved through their approved symbol.
        self.n_symbol_fallback: int = 0

    def parse(self) -> Dict[str, Set[str]]:
        gene_terms, symbols = parse_syngo_annotations(self.zip_path)
        index = parse_gene_accession_index(self.dat_path)
        gene_map, self.n_symbol_fallback = resolve_hgnc_accessions(
            gene_terms, symbols, index
        )
        remapped, self.coverage = remap_gene_annotations(
            gene_terms, gene_map, label="HGNC→UniProt (SynGO)"
        )
        return remapped
