"""Gene id → UniProt accession mapping, for gene-keyed annotation layers.

The statistics engine joins domains and terms on UniProt accessions
(``protein2ipr``'s key space), but some annotation databases are keyed by
*gene*: HPO's ``genes_to_phenotype.txt`` by NCBI GeneID, SynGO's annotations by
HGNC id. Those layers must re-key gene → accession at parse time, exactly as
the DOID layer re-keys OMIM → DOID — only on the protein axis of the map
instead of the term axis.

The translations already sit in the Swiss-Prot flat file the UniProt-native
layers read (``DR   GeneID; 1017; -.``, ``DR   HGNC; HGNC:1771; CDK2.``), so no
extra idmapping download is needed: one pass over the flat file builds all
three indexes (GeneID, HGNC id, HGNC-approved gene symbol) at once.

Mapping policy, mirroring :mod:`src.disease_ontology`:

* **Unmapped gene** (no reviewed UniProt entry cross-references it) — dropped,
  counted and logged, never silently discarded.
* **One-to-many** (one gene id cross-referenced by several accessions — real
  for readthrough loci and unresolved paralogs) — kept as a genuine expansion:
  the term goes to *all* of them, since choosing one arbitrarily is not
  reproducible.
* Coverage accounting *reuses* the generic :func:`src.remap.remap_values` on
  the inverted map (terms become the keys, gene ids the values being remapped)
  rather than duplicating its audited counting loop. The returned
  :class:`~src.remap.RemapCoverage` is axis-neutral, so its fields stay
  truthful here: the "value" counters range over gene ids, the "key" counters
  over ontology terms.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Set, Tuple

from loguru import logger

from src.remap import RemapCoverage, remap_values
from src.uniprot_annotation_source import iter_uniprot_entries


@dataclass(frozen=True)
class GeneAccessionMap:
    """One gene-id space → UniProt accessions, with its audit counts.

    Attributes:
        id_space: which ids the keys are (``"GeneID"``, ``"HGNC"``,
            ``"gene symbol"``).
        source_to_accessions: the mapping itself.
        n_entries: flat-file entries scanned to build it.
    """

    id_space: str
    source_to_accessions: Dict[str, Set[str]]
    n_entries: int = 0

    def __len__(self) -> int:
        return len(self.source_to_accessions)

    def targets(self, source_id: str) -> Set[str]:
        """Accessions cross-referencing ``source_id`` (empty if unmapped)."""
        return self.source_to_accessions.get(source_id, set())

    @property
    def n_one_to_many(self) -> int:
        """Gene ids carried by more than one accession."""
        return sum(1 for accs in self.source_to_accessions.values() if len(accs) > 1)


@dataclass(frozen=True)
class GeneAccessionIndex:
    """Every gene-id index one flat-file pass yields.

    Attributes:
        geneid: NCBI GeneID → accessions (``DR   GeneID; 1017; -.``).
        hgnc: HGNC id → accessions (``DR   HGNC; HGNC:1771; CDK2.``).
        symbol: HGNC-approved symbol → accessions (the third field of the same
            ``DR HGNC`` line), the fallback for annotations that only carry a
            symbol.
    """

    geneid: GeneAccessionMap
    hgnc: GeneAccessionMap
    symbol: GeneAccessionMap


def parse_gene_accession_index(dat_path: Path) -> GeneAccessionIndex:
    """Build every gene→accession index from one pass over the flat file.

    Gene ids are species-scoped by construction (NCBI GeneIDs are unique across
    species; HGNC is human-only), so no explicit species filter is needed: a
    human gene id can only match the human entry that cross-references it.
    """
    logger.info(f"Building gene → accession indexes from {dat_path}")
    geneid: Dict[str, Set[str]] = defaultdict(set)
    hgnc: Dict[str, Set[str]] = defaultdict(set)
    symbol: Dict[str, Set[str]] = defaultdict(set)
    n_entries = 0
    for entry in iter_uniprot_entries(Path(dat_path)):
        n_entries += 1
        if entry.accession is None:
            continue
        for db, xref_id, xref_type in entry.cross_refs:
            # "-" is UniProt's placeholder field; skipping it mirrors the
            # term filter parse_uniprot_cross_refs applies (ids never carry it
            # in practice, but the raw iterator leaves cleaning to us).
            if xref_id == "-":
                continue
            if db == "GeneID":
                geneid[xref_id].add(entry.accession)
            elif db == "HGNC":
                hgnc[xref_id].add(entry.accession)
                if xref_type and xref_type != "-":
                    symbol[xref_type].add(entry.accession)
    index = GeneAccessionIndex(
        geneid=GeneAccessionMap("GeneID", dict(geneid), n_entries),
        hgnc=GeneAccessionMap("HGNC", dict(hgnc), n_entries),
        symbol=GeneAccessionMap("gene symbol", dict(symbol), n_entries),
    )
    logger.info(
        f"  Entries scanned: {n_entries:,}; GeneID ids: {len(index.geneid):,} "
        f"({index.geneid.n_one_to_many:,} one-to-many); HGNC ids: "
        f"{len(index.hgnc):,}; symbols: {len(index.symbol):,}"
    )
    return index


def _invert(mapping: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    inverted: Dict[str, Set[str]] = defaultdict(set)
    for key, values in mapping.items():
        for value in values:
            inverted[value].add(key)
    return dict(inverted)


def remap_gene_annotations(
    gene_terms: Dict[str, Set[str]],
    gene_map: GeneAccessionMap,
    label: str,
) -> Tuple[Dict[str, Set[str]], RemapCoverage]:
    """Re-key a ``{gene id: {term}}`` map onto ``{accession: {term}}``.

    Implemented as the generic :func:`src.remap.remap_values` on the inverted
    map, so the unmapped/one-to-many accounting is the same audited code the
    DOID layer uses. The coverage fields read naturally: its *values* are the
    gene ids being remapped (``value_coverage`` is the fraction of genes that
    mapped, ``unmapped_values`` the genes that mapped to nothing) and its
    *keys* are the ontology terms.
    """
    remapped, coverage = remap_values(
        _invert(gene_terms),
        gene_map,
        label,
        key_label="term",
        value_label="gene",
        target_label="reviewed UniProt accession",
    )
    return _invert(remapped), coverage
