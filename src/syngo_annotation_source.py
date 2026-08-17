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
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

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


def _iter_sheet_rows(
    zip_path: Path, member: str, required: Iterable[str]
) -> Iterator[Dict[str, object]]:
    """Yield each data row of one xlsx inside the SynGO zip, keyed by header.

    The zip does not name a sheet, so the first sheet whose header row carries
    every ``required`` column is parsed; a release that renames a column (or
    ships an empty sheet) raises :class:`ValueError` naming the member and the
    missing columns instead of yielding a silently empty layer.

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
    try:
        seen_headers: Dict[str, List[str]] = {}
        for sheet in workbook.worksheets:
            rows = sheet.iter_rows(values_only=True)
            # next(rows, None), not next(rows): inside a generator a bare
            # StopIteration from an empty sheet would surface as an opaque
            # RuntimeError (PEP 479) rather than the ValueError below.
            first = next(rows, None)
            if first is None:
                continue
            header = [str(cell) for cell in first]
            seen_headers[sheet.title] = header
            if set(required) <= set(header):
                for row in rows:
                    yield dict(zip(header, row))
                return
        raise ValueError(
            f"{member} in {zip_path}: no sheet carries the expected "
            f"column(s) {sorted(required)}; "
            + (
                f"headers found: {seen_headers}"
                if seen_headers
                else "every sheet is empty"
            )
        )
    finally:
        workbook.close()


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
    for row in _iter_sheet_rows(
        zip_path, _ANNOTATIONS_MEMBER, required=("hgnc_id", "hgnc_symbol", "go_id")
    ):
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
    for row in _iter_sheet_rows(
        zip_path, _ONTOLOGIES_MEMBER, required=("id", "parent_id")
    ):
        term = _cell(row, "id")
        parent = _cell(row, "parent_id")
        if term and parent:
            child_to_parents.setdefault(term, set()).add(parent)
    return child_to_parents


@dataclass(frozen=True)
class HGNCResolution:
    """How the annotated HGNC ids resolved to accessions, with audit counts.

    Attributes:
        gene_map: the composite map, restricted to the annotated genes.
        n_genes: annotated HGNC ids considered.
        n_by_id: resolved by their HGNC id (``DR HGNC`` match).
        n_by_symbol: resolved by approved symbol, with no contradicting id.
        n_symbol_conflicts: the symbol matched an entry, but every matched
            entry cross-references a *different* HGNC id — a stale or
            reassigned symbol that would bind another gene's protein, so the
            match is rejected rather than trusted.
        n_unresolved: genes no route resolved (these are the ids that then
            appear in the remap coverage's ``unmapped_values``).
    """

    gene_map: GeneAccessionMap
    n_genes: int = 0
    n_by_id: int = 0
    n_by_symbol: int = 0
    n_symbol_conflicts: int = 0
    n_unresolved: int = 0


def resolve_hgnc_accessions(
    gene_terms: Dict[str, Set[str]],
    symbols: Dict[str, str],
    index: GeneAccessionIndex,
) -> HGNCResolution:
    """Resolve each annotated HGNC id to accessions, symbol as guarded fallback.

    A symbol-only match is cross-checked against the matched entry's own
    ``DR HGNC`` id: if the entry names a different gene, the symbol is stale
    (an HGNC merge) or reassigned (another gene now owns it), and without HGNC
    history the two cannot be told apart — so the match is rejected and
    counted (``n_symbol_conflicts``) rather than silently binding a different
    gene's protein. Since the symbol index is itself built from ``DR HGNC``
    lines, this guard makes the fallback deliberately conservative.
    """
    accession_ids: Dict[str, Set[str]] = defaultdict(set)
    for hgnc_id, accessions in index.hgnc.source_to_accessions.items():
        for accession in accessions:
            accession_ids[accession].add(hgnc_id)

    resolved: Dict[str, Set[str]] = {}
    n_by_id = 0
    n_by_symbol = 0
    n_conflicts = 0
    for gene in gene_terms:
        accessions = index.hgnc.targets(gene)
        if accessions:
            n_by_id += 1
        else:
            symbol = symbols.get(gene, "")
            candidates = index.symbol.targets(symbol) if symbol else set()
            accessions = {
                accession
                for accession in candidates
                if not accession_ids.get(accession) or gene in accession_ids[accession]
            }
            if accessions:
                n_by_symbol += 1
            elif candidates:
                n_conflicts += 1
        if accessions:
            resolved[gene] = accessions

    resolution = HGNCResolution(
        gene_map=GeneAccessionMap("HGNC(+symbol)", resolved, index.hgnc.n_entries),
        n_genes=len(gene_terms),
        n_by_id=n_by_id,
        n_by_symbol=n_by_symbol,
        n_symbol_conflicts=n_conflicts,
        n_unresolved=len(gene_terms) - n_by_id - n_by_symbol,
    )
    logger.info(
        f"  HGNC resolution: {resolution.n_by_id:,} by id, "
        f"{resolution.n_by_symbol:,} by symbol, "
        f"{resolution.n_symbol_conflicts:,} symbol conflicts rejected, "
        f"{resolution.n_unresolved:,} unresolved of {resolution.n_genes:,} genes"
    )
    return resolution


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
        #: Populated by :meth:`parse`: how each HGNC id resolved (by id, by
        #: guarded symbol fallback, rejected symbol conflict, or not at all).
        self.resolution: Optional[HGNCResolution] = None

    def parse(self) -> Dict[str, Set[str]]:
        gene_terms, symbols = parse_syngo_annotations(self.zip_path)
        index = parse_gene_accession_index(self.dat_path)
        self.resolution = resolve_hgnc_accessions(gene_terms, symbols, index)
        remapped, self.coverage = remap_gene_annotations(
            gene_terms, self.resolution.gene_map, label="HGNC→UniProt (SynGO)"
        )
        return remapped
