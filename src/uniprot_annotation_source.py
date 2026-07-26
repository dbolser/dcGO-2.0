"""UniProt-native annotation sources: harvest terms that UniProt already keys by accession.

UniProt accessions are the protein universe here — ``protein2ipr`` domains, GOA,
and Expasy ENZYME are all keyed by them. So the cheapest ontologies to add are
the ones UniProt *already* carries per accession, needing no identifier mapping:

* **DR (database cross-reference) lines** point each entry at external resources
  — ``Reactome``, ``KEGG``, ``GO``, ``InterPro``, ``MIM``, ``Orphanet``,
  ``DisGeNET``, ``DrugBank``, ``ChEMBL``, ``PANTHER`` … One parser, many
  ontologies: pick the database name.
* **KW (keyword) lines** are a controlled vocabulary (the UniProt keyword list)
  spanning function, disease, biological process, and more.

This module parses the UniProt flat file (``uniprot_sprot.dat.gz``) and exposes
those as :class:`AnnotationSource` implementations, so the dcGO engine can
associate domains with any UniProt-native vocabulary the same way it does GO.

Flat-file entry shape (``//``-delimited)::

    AC   P07327; B2R5V5;
    DR   Reactome; R-HSA-71384; Ethanol oxidation.
    DR   KEGG; hsa:124; .
    KW   Metal-binding; NAD; Oxidoreductase; Zinc.
    //
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Set, Tuple

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec

# Reactome stable ids look like "R-HSA-71384"; keywords carry no id prefix.
REACTOME_SPEC = OntologySpec(
    ontology_id="Reactome", name="Reactome pathways", term_prefix="R-"
)
KEYWORD_SPEC = OntologySpec(ontology_id="UniProtKW", name="UniProt keywords")
# OMIM MIM numbers (numeric, no prefix). UniProt DR MIM lines are typed
# gene/phenotype; disease association uses the phenotype entries.
DISEASE_SPEC = OntologySpec(ontology_id="OMIM", name="OMIM disease (phenotype)")


def _open_text(path: Path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"UniProt flat file not found: {path}")
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def _split_keywords(kw_parts: List[str]) -> List[str]:
    """Join wrapped ``KW`` line payloads and split into individual keywords.

    UniProt only ever wraps ``KW`` lines after a ``;``, so joining and splitting
    on ``;`` is safe. The final keyword ends with ``.``.
    """
    if not kw_parts:
        return []
    text = " ".join(p.strip() for p in kw_parts).rstrip().rstrip(".")
    return [k.strip() for k in text.split(";") if k.strip()]


def _iter_entries(
    path: Path,
) -> Iterator[Tuple[str | None, List[Tuple[str, str, str]], List[str]]]:
    """Yield ``(primary_accession, dr_triples, keywords)`` per flat-file entry.

    ``dr_triples`` is a list of ``(database, external_id, id_type)``, where
    ``id_type`` is the DR line's third field (e.g. ``"gene"``/``"phenotype"`` for
    ``MIM``, ``""`` when absent). ``keywords`` is the entry's keyword list. The
    primary accession is the first accession on the first ``AC`` line — the same
    key space as ``protein2ipr``.
    """
    with _open_text(path) as f:
        accession: str | None = None
        dr_triples: List[Tuple[str, str, str]] = []
        kw_parts: List[str] = []

        for line in f:
            tag = line[:2]
            if tag == "AC":
                if accession is None:
                    first = line[5:].split(";")[0].strip()
                    accession = first or None
            elif tag == "DR":
                fields = line[5:].split(";")
                if len(fields) >= 2:
                    db = fields[0].strip()
                    xref_id = fields[1].strip()
                    id_type = fields[2].strip().rstrip(".") if len(fields) >= 3 else ""
                    if db and xref_id:
                        dr_triples.append((db, xref_id, id_type))
            elif tag == "KW":
                kw_parts.append(line[5:])
            elif line.startswith("//"):
                yield accession, dr_triples, _split_keywords(kw_parts)
                accession = None
                dr_triples = []
                kw_parts = []


def parse_uniprot_cross_refs(
    path: Path, database: str, id_type: str | None = None
) -> Dict[str, Set[str]]:
    """Return ``{accession: {external_id}}`` for one DR database (e.g. ``"Reactome"``).

    Args:
        path: UniProt flat file.
        database: exact DR database name (``"Reactome"``, ``"MIM"``, …).
        id_type: if given, keep only cross-references whose third DR field
            matches (e.g. ``"phenotype"`` to select OMIM disease entries and
            drop the ``"gene"`` ones). ``None`` keeps all.
    """
    label = database if id_type is None else f"{database}/{id_type}"
    logger.info(f"Parsing UniProt cross-references ({label}) from {path}")
    result: Dict[str, Set[str]] = defaultdict(set)
    n_entries = 0
    for accession, dr_triples, _kw in _iter_entries(path):
        n_entries += 1
        if accession is None:
            continue
        for db, xref_id, xref_type in dr_triples:
            if db == database and (id_type is None or xref_type == id_type):
                result[accession].add(xref_id)
    logger.info(
        f"  Entries scanned: {n_entries:,}; proteins with {label}: {len(result):,}"
    )
    return dict(result)


def parse_uniprot_keywords(path: Path) -> Dict[str, Set[str]]:
    """Return ``{accession: {keyword}}`` from UniProt ``KW`` lines."""
    logger.info(f"Parsing UniProt keywords from {path}")
    result: Dict[str, Set[str]] = defaultdict(set)
    n_entries = 0
    for accession, _dr, keywords in _iter_entries(path):
        n_entries += 1
        if accession is None or not keywords:
            continue
        result[accession].update(keywords)
    logger.info(
        f"  Entries scanned: {n_entries:,}; proteins with keywords: {len(result):,}"
    )
    return dict(result)


class UniProtCrossRefAnnotationSource(AnnotationSource):
    """Domain annotations from one UniProt DR cross-reference database.

    ``database`` is the exact DR database name as it appears in the flat file
    (``"Reactome"``, ``"KEGG"``, ``"MIM"``, …). Because the flat file is keyed by
    UniProt accession, the resulting terms join directly to the domain data.
    """

    def __init__(
        self,
        dat_path: Path,
        database: str,
        spec: OntologySpec,
        id_type: str | None = None,
    ) -> None:
        self.dat_path = Path(dat_path)
        self.database = database
        self.spec = spec
        self.id_type = id_type

    def parse(self) -> Dict[str, Set[str]]:
        return parse_uniprot_cross_refs(
            self.dat_path, self.database, id_type=self.id_type
        )


class UniProtKeywordAnnotationSource(AnnotationSource):
    """Domain annotations from UniProt keywords (``KW`` lines)."""

    def __init__(self, dat_path: Path, spec: OntologySpec = KEYWORD_SPEC) -> None:
        self.dat_path = Path(dat_path)
        self.spec = spec

    def parse(self) -> Dict[str, Set[str]]:
        return parse_uniprot_keywords(self.dat_path)


def reactome_source(dat_path: Path) -> UniProtCrossRefAnnotationSource:
    """Convenience factory for a Reactome-pathway annotation source."""
    return UniProtCrossRefAnnotationSource(dat_path, "Reactome", REACTOME_SPEC)


def disease_source(dat_path: Path) -> UniProtCrossRefAnnotationSource:
    """Convenience factory for an OMIM disease annotation source.

    Uses UniProt ``DR MIM`` cross-references restricted to ``phenotype`` entries
    (dropping the ``gene`` MIM links), i.e. the disease side of OMIM.
    """
    return UniProtCrossRefAnnotationSource(
        dat_path, "MIM", DISEASE_SPEC, id_type="phenotype"
    )


def parse_reactome_relations(
    path: Path, species_prefix: str | None = None
) -> Dict[str, Set[str]]:
    """Parse Reactome ``ReactomePathwaysRelation.txt`` into ``{child: {parents}}``.

    The file is two tab-separated columns, ``parent_id<TAB>child_id`` (stable ids
    like ``R-HSA-71384``). Feed the result to
    :func:`src.hierarchy.closure_ancestors` to propagate domain→pathway
    associations up the pathway hierarchy.

    Args:
        path: the relations file (optionally gzipped).
        species_prefix: if given (e.g. ``"R-HSA-"``), keep only edges whose ids
            both start with it. ``None`` keeps all species.
    """
    child_to_parents: Dict[str, Set[str]] = defaultdict(set)
    with _open_text(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            parent, child = parts[0].strip(), parts[1].strip()
            if not parent or not child:
                continue
            if species_prefix and not (
                parent.startswith(species_prefix) and child.startswith(species_prefix)
            ):
                continue
            child_to_parents[child].add(parent)
    logger.info(
        f"Parsed Reactome hierarchy: {len(child_to_parents):,} pathways with parents"
    )
    return dict(child_to_parents)


def parse_keyword_hierarchy(path: Path) -> Dict[str, Set[str]]:
    """Parse the UniProt keyword list (``keywlist.txt``) into ``{keyword: {parents}}``.

    Keyword names are the terms our ``KW`` harvesting produces, and the keyword
    list encodes the hierarchy on ``HI`` lines::

        ID   2Fe-2S.
        HI   Ligand: Iron; Iron-sulfur; 2Fe-2S.

    Each ``HI`` line is a path ``Category: parent; …; thisKeyword``; the term
    immediately before the current keyword is its (a) direct parent. Keywords
    form a DAG (multiple ``HI`` lines → multiple parents). Feed the result to
    :func:`src.hierarchy.closure_ancestors`.
    """
    child_to_parents: Dict[str, Set[str]] = defaultdict(set)
    current: str | None = None
    with _open_text(path) as f:
        for line in f:
            tag = line[:2]
            if tag == "ID":
                current = line[5:].strip().rstrip(".") or None
            elif tag == "HI" and current is not None:
                # Drop the "Category:" prefix, then split the path on ';'.
                payload = line[5:].split(":", 1)
                path_part = payload[1] if len(payload) > 1 else payload[0]
                items = [
                    x.strip().rstrip(".") for x in path_part.split(";") if x.strip()
                ]
                # Path ends at the current keyword; its parent is the one before.
                if len(items) >= 2 and items[-1] == current:
                    child_to_parents[current].add(items[-2])
            elif line.startswith("//"):
                current = None
    logger.info(
        f"Parsed keyword hierarchy: {len(child_to_parents):,} keywords with parents"
    )
    return dict(child_to_parents)
