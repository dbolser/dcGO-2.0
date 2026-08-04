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
* **CC (comment) blocks** carry curated annotation against further controlled
  vocabularies: ``SUBCELLULAR LOCATION`` (the ``subcell.txt`` ``SL-`` terms),
  ``COFACTOR`` (ChEBI) and ``CATALYTIC ACTIVITY`` (Rhea reactions, a finer
  enzymology layer than EC).
* **FT (feature) qualifiers** name the chemistry bound at annotated residues —
  ``/ligand_id="ChEBI:…"`` — which is the layer closest to what a domain is.

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
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple

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
# Annotation layers that are *not* DR cross-references: they are curated into
# the entry body (CC comment blocks and FT feature qualifiers) against
# controlled vocabularies of their own.
SUBCELLULAR_SPEC = OntologySpec(
    ontology_id="SL", name="UniProt subcellular location", term_prefix="SL-"
)
LIGAND_SPEC = OntologySpec(
    ontology_id="ChEBI-ligand", name="Bound ligand (ChEBI)", term_prefix="CHEBI:"
)
COFACTOR_SPEC = OntologySpec(
    ontology_id="ChEBI-cofactor", name="Cofactor (ChEBI)", term_prefix="CHEBI:"
)
RHEA_SPEC = OntologySpec(
    ontology_id="Rhea", name="Catalysed reaction (Rhea)", term_prefix="RHEA:"
)

# "{ECO:0000269|PubMed:2738060}" evidence tags, stripped before parsing prose.
_EVIDENCE_RE = re.compile(r"\{[^{}]*\}")
# FT qualifier: /ligand_id="ChEBI:CHEBI:29105"
_LIGAND_ID_RE = re.compile(r'/ligand_id="ChEBI:(CHEBI:\d+)"')
# CC COFACTOR / CATALYTIC ACTIVITY cross-references.
_CHEBI_RE = re.compile(r"ChEBI:(CHEBI:\d+)")
_RHEA_RE = re.compile(r"Rhea:(RHEA:\d+)")


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


@dataclass
class UniProtEntry:
    """One flat-file entry, reduced to the annotation layers dcGO can use.

    Attributes:
        accession: primary accession (first id on the first ``AC`` line) — the
            same key space as ``protein2ipr``.
        cross_refs: ``(database, external_id, id_type)`` per ``DR`` line, where
            ``id_type`` is the third field (``"gene"``/``"phenotype"`` for
            ``MIM``, ``""`` when absent).
        keywords: the entry's ``KW`` controlled-vocabulary terms.
        cc_blocks: ``(topic, payload)`` per ``CC   -!- TOPIC: …`` comment block,
            continuation lines joined with spaces. Only collected when
            ``want_cc`` is set.
        ligand_ids: ChEBI ids from ``FT`` ``/ligand_id`` qualifiers (bound
            ligands, cofactors and substrates at annotated sites). Only
            collected when ``want_ft`` is set.
    """

    accession: Optional[str] = None
    cross_refs: List[Tuple[str, str, str]] = field(default_factory=list)
    keywords: List[str] = field(default_factory=list)
    cc_blocks: List[Tuple[str, str]] = field(default_factory=list)
    ligand_ids: Set[str] = field(default_factory=set)


def _iter_entries(
    path: Path, want_cc: bool = False, want_ft: bool = False
) -> Iterator[UniProtEntry]:
    """Yield one :class:`UniProtEntry` per flat-file entry.

    ``want_cc``/``want_ft`` are opt-in because the flat file is ~1 GB
    compressed: a Reactome or keyword run should not pay to assemble comment
    blocks it will never read.
    """
    with _open_text(path) as f:
        entry = UniProtEntry()
        cc_topic: str | None = None
        cc_parts: List[str] = []

        def close_cc_block() -> None:
            if cc_topic is not None:
                entry.cc_blocks.append((cc_topic, " ".join(cc_parts).strip()))

        for line in f:
            tag = line[:2]
            if tag == "AC":
                if entry.accession is None:
                    first = line[5:].split(";")[0].strip()
                    entry.accession = first or None
            elif tag == "DR":
                fields = line[5:].split(";")
                if len(fields) >= 2:
                    db = fields[0].strip()
                    xref_id = fields[1].strip()
                    id_type = fields[2].strip().rstrip(".") if len(fields) >= 3 else ""
                    if db and xref_id:
                        entry.cross_refs.append((db, xref_id, id_type))
            elif tag == "KW":
                entry.keywords.append(line[5:])
            elif tag == "CC" and want_cc:
                payload = line[5:].rstrip("\n")
                if payload.startswith("-!- "):
                    close_cc_block()
                    topic, _, rest = payload[4:].partition(":")
                    cc_topic, cc_parts = topic.strip(), [rest.strip()]
                elif payload.startswith("---"):
                    # The copyright block every entry ends with is CC-tagged but
                    # is not a comment topic; without this it would be swallowed
                    # as continuation text by whichever topic came last.
                    close_cc_block()
                    cc_topic, cc_parts = None, []
                elif cc_topic is not None:
                    cc_parts.append(payload.strip())
            elif tag == "FT" and want_ft:
                match = _LIGAND_ID_RE.search(line)
                if match:
                    entry.ligand_ids.add(match.group(1))
            elif line.startswith("//"):
                close_cc_block()
                entry.keywords = _split_keywords(entry.keywords)
                yield entry
                entry = UniProtEntry()
                cc_topic, cc_parts = None, []


def parse_uniprot_accessions(path: Path) -> Set[str]:
    """Every primary accession in a UniProt flat file.

    Temporal evaluations need to know which proteins *existed* in a snapshot,
    not merely which carried a given annotation: a protein absent from the t0
    release cannot have been "predicted" to gain a term, it simply had no entry
    yet. An :class:`AnnotationSource` cannot answer that, because it only
    reports proteins that carry terms of its own ontology.
    """
    logger.info(f"Reading the accession set from {path}")
    accessions: Set[str] = set()
    for entry in _iter_entries(path):
        if entry.accession:
            accessions.add(entry.accession)
    logger.info(f"  Entries: {len(accessions):,}")
    return accessions


def parse_uniprot_cross_refs(
    path: Path,
    database: str,
    id_type: str | None = None,
    term_from_id_type: bool = False,
) -> Dict[str, Set[str]]:
    """Return ``{accession: {external_id}}`` for one DR database (e.g. ``"Reactome"``).

    Args:
        path: UniProt flat file.
        database: exact DR database name (``"Reactome"``, ``"MIM"``, …).
        id_type: if given, keep only cross-references whose third DR field
            matches (e.g. ``"phenotype"`` to select OMIM disease entries and
            drop the ``"gene"`` ones). ``None`` keeps all.
        term_from_id_type: take the *third* DR field as the term instead of the
            id. A handful of databases key the DR line by the protein itself and
            carry the classification in that field — Pharos target development
            level (``DR Pharos; P31946; Tbio.``), CD-CODE condensate name — so
            the vocabulary lives there, not in the id.
    """
    label = database if id_type is None else f"{database}/{id_type}"
    logger.info(f"Parsing UniProt cross-references ({label}) from {path}")
    result: Dict[str, Set[str]] = defaultdict(set)
    n_entries = 0
    for entry in _iter_entries(path):
        n_entries += 1
        if entry.accession is None:
            continue
        for db, xref_id, xref_type in entry.cross_refs:
            if db == database and (id_type is None or xref_type == id_type):
                term = xref_type if term_from_id_type else xref_id
                if term and term != "-":
                    result[entry.accession].add(term)
    logger.info(
        f"  Entries scanned: {n_entries:,}; proteins with {label}: {len(result):,}"
    )
    return dict(result)


def parse_uniprot_keywords(path: Path) -> Dict[str, Set[str]]:
    """Return ``{accession: {keyword}}`` from UniProt ``KW`` lines."""
    logger.info(f"Parsing UniProt keywords from {path}")
    result: Dict[str, Set[str]] = defaultdict(set)
    n_entries = 0
    for entry in _iter_entries(path):
        n_entries += 1
        if entry.accession is None or not entry.keywords:
            continue
        result[entry.accession].update(entry.keywords)
    logger.info(
        f"  Entries scanned: {n_entries:,}; proteins with keywords: {len(result):,}"
    )
    return dict(result)


@dataclass
class SubcellVocabulary:
    """The UniProt subcellular-location controlled vocabulary (``subcell.txt``).

    Attributes:
        accession_of: normalised location string → ``SL-`` accession. Keys cover
            the ``SL`` content line (what ``CC   -!- SUBCELLULAR LOCATION``
            lines actually contain), the entry name (``ID``/``IT``/``IO``) and
            every ``SY`` synonym.
        name_of: ``SL-`` accession → canonical entry name.
        child_to_parents: ``SL-`` accession → direct parents, from the ``HI``
            (is-a) and ``HP`` (part-of) lines — both are traversed by the True
            Path Rule.
    """

    accession_of: Dict[str, str] = field(default_factory=dict)
    name_of: Dict[str, str] = field(default_factory=dict)
    child_to_parents: Dict[str, Set[str]] = field(default_factory=dict)


def _normalise_location(text: str) -> str:
    """Canonical form for matching a location string against the vocabulary."""
    return text.strip().rstrip(".").strip().casefold()


def parse_subcell_vocabulary(path: Path) -> SubcellVocabulary:
    """Parse UniProt ``subcell.txt`` into a :class:`SubcellVocabulary`.

    Entry shape (locations use ``ID``, membrane topologies ``IT``, orientations
    ``IO``; all carry an ``AC`` and an ``SL`` content line)::

        ID   Nucleolus.
        AC   SL-0188
        SY   Nucleoli.
        SL   Nucleus, nucleolus.
        HP   Nucleus.
        //
    """
    vocab = SubcellVocabulary()
    parents_by_name: Dict[str, Set[str]] = defaultdict(set)

    name: str | None = None
    accession: str | None = None
    aliases: List[str] = []
    parent_names: List[str] = []

    with _open_text(path) as f:
        for line in f:
            tag = line[:2]
            payload = line[5:].strip()
            if tag in ("ID", "IT", "IO"):
                name = payload.rstrip(".")
                aliases = [name]
            elif tag == "AC":
                accession = payload
            elif tag == "SY" and name is not None:
                aliases.extend(s.strip().rstrip(".") for s in payload.split(";"))
            elif tag == "SL" and name is not None:
                aliases.append(payload.rstrip("."))
            elif tag in ("HI", "HP") and name is not None:
                parent_names.append(payload.rstrip("."))
            elif line.startswith("//"):
                if name and accession:
                    vocab.name_of[accession] = name
                    for alias in aliases:
                        if alias:
                            vocab.accession_of.setdefault(
                                _normalise_location(alias), accession
                            )
                    if parent_names:
                        parents_by_name[accession] = set(parent_names)
                name, accession, aliases, parent_names = None, None, [], []

    # HI/HP reference parents by *name*; re-key them to accessions.
    for child_ac, names in parents_by_name.items():
        parents = {
            vocab.accession_of[key]
            for key in (_normalise_location(n) for n in names)
            if key in vocab.accession_of
        }
        if parents:
            vocab.child_to_parents[child_ac] = parents

    logger.info(
        f"Parsed subcellular-location vocabulary: {len(vocab.name_of):,} terms, "
        f"{len(vocab.child_to_parents):,} with parents"
    )
    return vocab


def _split_location_statements(payload: str) -> List[str]:
    """Split a ``SUBCELLULAR LOCATION`` comment into individual location strings.

    The comment is a run of ``.``-terminated statements, each optionally
    prefixed with an isoform tag and optionally carrying ``;``-separated
    topology/orientation qualifiers, with a free-text ``Note=`` tail::

        [Isoform 1]: Cell membrane; Single-pass type II membrane protein. Nucleus. Note=…

    → ``["Cell membrane", "Single-pass type II membrane protein", "Nucleus"]``
    """
    text = _EVIDENCE_RE.sub("", payload).split("Note=")[0]
    pieces: List[str] = []
    for statement in text.split("."):
        statement = statement.strip()
        if not statement:
            continue
        if statement.startswith("["):  # "[Isoform 1]: Nucleus"
            statement = statement.partition("]:")[2].strip()
        for piece in statement.split(";"):
            piece = piece.strip()
            if piece:
                pieces.append(piece)
    return pieces


def parse_subcellular_locations(path: Path, subcell_path: Path) -> Dict[str, Set[str]]:
    """Return ``{accession: {SL-id}}`` from ``CC   -!- SUBCELLULAR LOCATION`` blocks.

    Locations are curated free-standing prose *against a controlled vocabulary*:
    each statement matches an entry in ``subcell.txt``, so mapping them to
    ``SL-`` accessions turns the comment block into a proper ontology layer
    (with a hierarchy, via :func:`parse_subcell_vocabulary`).

    Args:
        path: UniProt flat file.
        subcell_path: ``subcell.txt`` controlled vocabulary.
    """
    logger.info(f"Parsing UniProt subcellular locations from {path}")
    vocab = parse_subcell_vocabulary(subcell_path)

    result: Dict[str, Set[str]] = defaultdict(set)
    unmatched: Dict[str, int] = defaultdict(int)
    n_entries = 0
    for entry in _iter_entries(path, want_cc=True):
        n_entries += 1
        if entry.accession is None:
            continue
        for topic, payload in entry.cc_blocks:
            if topic != "SUBCELLULAR LOCATION":
                continue
            for piece in _split_location_statements(payload):
                accession = vocab.accession_of.get(_normalise_location(piece))
                if accession:
                    result[entry.accession].add(accession)
                else:
                    unmatched[piece] += 1

    logger.info(
        f"  Entries scanned: {n_entries:,}; proteins with a location: {len(result):,}"
    )
    if unmatched:
        worst = sorted(unmatched.items(), key=lambda kv: -kv[1])[:3]
        logger.warning(
            f"  {sum(unmatched.values()):,} location strings did not match "
            f"subcell.txt ({len(unmatched):,} distinct); most frequent: "
            + ", ".join(f"{text!r} ×{n}" for text, n in worst)
        )
    return dict(result)


def parse_binding_ligands(path: Path) -> Dict[str, Set[str]]:
    """Return ``{accession: {ChEBI id}}`` from ``FT`` ``/ligand_id`` qualifiers.

    UniProt annotates the chemistry a protein binds at the residue level
    (``FT   BINDING`` and friends), naming each ligand with a ChEBI id. Binding
    sites are exactly the kind of feature a *domain* carries, so this is the
    most domain-proximal of the UniProt-native layers.
    """
    logger.info(f"Parsing UniProt binding-site ligands (ChEBI) from {path}")
    result: Dict[str, Set[str]] = defaultdict(set)
    n_entries = 0
    for entry in _iter_entries(path, want_ft=True):
        n_entries += 1
        if entry.accession is None or not entry.ligand_ids:
            continue
        result[entry.accession] |= entry.ligand_ids
    logger.info(
        f"  Entries scanned: {n_entries:,}; proteins with a ligand: {len(result):,}"
    )
    return dict(result)


def _parse_cc_xrefs(
    path: Path, topic: str, pattern: re.Pattern, label: str
) -> Dict[str, Set[str]]:
    """Return ``{accession: {id}}`` for every ``pattern`` match in one CC topic."""
    logger.info(f"Parsing UniProt {label} from {path}")
    result: Dict[str, Set[str]] = defaultdict(set)
    n_entries = 0
    for entry in _iter_entries(path, want_cc=True):
        n_entries += 1
        if entry.accession is None:
            continue
        for block_topic, payload in entry.cc_blocks:
            if block_topic == topic:
                result[entry.accession].update(pattern.findall(payload))
    result = {acc: terms for acc, terms in result.items() if terms}
    logger.info(
        f"  Entries scanned: {n_entries:,}; proteins with {label}: {len(result):,}"
    )
    return result


def parse_cofactors(path: Path) -> Dict[str, Set[str]]:
    """Return ``{accession: {ChEBI id}}`` from ``CC   -!- COFACTOR`` blocks."""
    return _parse_cc_xrefs(path, "COFACTOR", _CHEBI_RE, "cofactors (ChEBI)")


def parse_catalysed_reactions(path: Path) -> Dict[str, Set[str]]:
    """Return ``{accession: {Rhea id}}`` from ``CC   -!- CATALYTIC ACTIVITY`` blocks.

    Rhea reaction ids are a finer-grained enzymology layer than EC: one EC
    number can cover many Rhea reactions differing in substrate.
    """
    return _parse_cc_xrefs(
        path, "CATALYTIC ACTIVITY", _RHEA_RE, "catalysed reactions (Rhea)"
    )


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
        term_from_id_type: bool = False,
    ) -> None:
        self.dat_path = Path(dat_path)
        self.database = database
        self.spec = spec
        self.id_type = id_type
        self.term_from_id_type = term_from_id_type

    def parse(self) -> Dict[str, Set[str]]:
        return parse_uniprot_cross_refs(
            self.dat_path,
            self.database,
            id_type=self.id_type,
            term_from_id_type=self.term_from_id_type,
        )


class UniProtKeywordAnnotationSource(AnnotationSource):
    """Domain annotations from UniProt keywords (``KW`` lines)."""

    def __init__(self, dat_path: Path, spec: OntologySpec = KEYWORD_SPEC) -> None:
        self.dat_path = Path(dat_path)
        self.spec = spec

    def parse(self) -> Dict[str, Set[str]]:
        return parse_uniprot_keywords(self.dat_path)


class UniProtSubcellularAnnotationSource(AnnotationSource):
    """Domain annotations from ``CC   -!- SUBCELLULAR LOCATION`` comments.

    Needs the ``subcell.txt`` controlled vocabulary to turn the curated location
    strings into stable ``SL-`` accessions (and to supply the hierarchy).
    """

    def __init__(
        self,
        dat_path: Path,
        subcell_path: Path,
        spec: OntologySpec = SUBCELLULAR_SPEC,
    ) -> None:
        self.dat_path = Path(dat_path)
        self.subcell_path = Path(subcell_path)
        self.spec = spec

    def parse(self) -> Dict[str, Set[str]]:
        return parse_subcellular_locations(self.dat_path, self.subcell_path)


class UniProtLigandAnnotationSource(AnnotationSource):
    """Domain annotations from ``FT`` ``/ligand_id`` ChEBI qualifiers."""

    def __init__(self, dat_path: Path, spec: OntologySpec = LIGAND_SPEC) -> None:
        self.dat_path = Path(dat_path)
        self.spec = spec

    def parse(self) -> Dict[str, Set[str]]:
        return parse_binding_ligands(self.dat_path)


class UniProtCofactorAnnotationSource(AnnotationSource):
    """Domain annotations from ``CC   -!- COFACTOR`` ChEBI cross-references."""

    def __init__(self, dat_path: Path, spec: OntologySpec = COFACTOR_SPEC) -> None:
        self.dat_path = Path(dat_path)
        self.spec = spec

    def parse(self) -> Dict[str, Set[str]]:
        return parse_cofactors(self.dat_path)


class UniProtRheaAnnotationSource(AnnotationSource):
    """Domain annotations from ``CC   -!- CATALYTIC ACTIVITY`` Rhea reactions."""

    def __init__(self, dat_path: Path, spec: OntologySpec = RHEA_SPEC) -> None:
        self.dat_path = Path(dat_path)
        self.spec = spec

    def parse(self) -> Dict[str, Set[str]]:
        return parse_catalysed_reactions(self.dat_path)


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
