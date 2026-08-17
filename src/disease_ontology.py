"""Disease Ontology: re-key UniProt's OMIM/Orphanet disease layer onto DOID terms.

``--ontology disease`` harvests UniProt ``DR MIM; …; phenotype.`` lines, so its
terms are raw OMIM phenotype ids. OMIM is a *catalogue*, not an ontology: the
ids carry no hierarchy at all, which costs the layer twice over.

1. **No True Path Rule.** ``--enable-true-path`` refuses, because there are no
   ancestors to propagate to.
2. **No pooling.** OMIM splits one disease across many allelic/locus-specific
   entries ("spinocerebellar ataxia 1/2/3/…", each its own MIM number), so the
   evidence for a domain is scattered over a long tail of terms that each fall
   below ``min_proteins_per_association``. Nothing reaches significance.

The Human Disease Ontology (DO) fixes both: it is a proper ``is_a`` DAG *and*
it cross-references OMIM (``xref: MIM:<id>``, note the ``MIM`` prefix, not
``OMIM``) and Orphanet (``xref: ORDO:<id>``). Translating each protein's OMIM
ids to DOID terms **at parse time** — before the contingency tables are built —
lets those sparse phenotypes pool into a better-supported DO class, and then
propagate up the DAG.

Parse time is the whole point. A post-hoc re-labelling of the output could only
rename terms that already reached significance; the pooling that gets them
there has to happen in the protein→term map that the Fisher engine consumes.

Mapping policy
--------------
The OMIM→DO mapping is *not* one-to-one, so every non-injective case is decided
here explicitly and counted at parse time (:class:`XrefMapping` for the table,
:class:`RemapCoverage` for what it did to the annotations):

* **Unmapped source id** (no DO term cross-references it) — dropped from the
  DOID layer, but counted and logged, never silently discarded. We cannot invent
  a DO class for it, and carrying the bare OMIM id through would put a term with
  no ancestors back into a layer whose entire purpose is the hierarchy. The
  ``disease`` ontology key keeps the raw OMIM layer available for comparison.
* **One-to-many** (one MIM id cross-referenced by several DO terms) — kept as a
  genuine one-to-many expansion: the protein gets *all* of them. DO splitting a
  MIM entry into two disease classes is a real curation statement, dropping
  either loses evidence, and choosing one arbitrarily is not reproducible. The
  DAG pools them again at their common ancestor. It is also rare: 11 of the
  6,145 numeric MIM ids in the 2026-07-31 release.
* **Obsolete DO term** — skipped, unless its ``replaced_by`` resolves
  (transitively, with a cycle guard) to a live term, which is then used.
  Obsolete stanzas are excluded from the ``is_a`` graph, so an obsolete target
  would be an orphan that can never propagate — the worst of both worlds.
  ``consider:`` is deliberately *not* followed: it is a suggestion for a curator,
  not an assertion of equivalence.
* **Non-numeric MIM xrefs** — DO also cross-references OMIM *phenotypic series*
  (``MIM:PS267700``, 332 of them). UniProt ``DR MIM`` lines only ever carry plain
  numeric ids, so these simply never match; they are counted separately rather
  than reported as failures.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, FrozenSet, Optional, Set, Tuple

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec
from src.remap import RemapCoverage, remap_values
from src.uniprot_annotation_source import parse_uniprot_cross_refs

#: Disease Ontology terms, e.g. ``DOID:9352``.
DOID_SPEC = OntologySpec(
    ontology_id="DOID", name="Human Disease Ontology", term_prefix="DOID:"
)

_ID_RE = re.compile(r"^id:\s*(\S+)")
_REPLACED_BY_RE = re.compile(r"^replaced_by:\s*(\S+)")
_XREF_RE = re.compile(r"^xref:\s*([^\s:]+):(\S+)")


@dataclass(frozen=True)
class XrefMapping:
    """A ``<source database> id → {DOID}`` translation table, with its audit trail.

    Attributes:
        xref_prefix: the OBO ``xref`` prefix this table was built from
            (``"MIM"`` for OMIM, ``"ORDO"`` for Orphanet).
        source_to_doids: the mapping itself. Targets are live (non-obsolete)
            DO terms only; obsolete ones have been resolved through
            ``replaced_by`` or dropped.
        n_terms: ``[Term]`` stanzas read.
        n_obsolete_terms: of those, how many were ``is_obsolete: true``.
        n_xrefs: ``xref: <prefix>:…`` lines seen (the raw cross-reference count).
        n_source_ids: distinct source ids seen on those lines, mapped or not.
        n_non_numeric: source ids that are not plain numbers — OMIM phenotypic
            series (``PS267700``), which UniProt ``DR`` lines never carry.
        n_one_to_many: source ids that map to more than one live DO term.
        n_obsolete_resolved: xrefs pointing at an obsolete term that
            ``replaced_by`` rescued.
        n_obsolete_dropped: xrefs pointing at an obsolete term with no usable
            replacement.
        source_ids_without_target: source ids seen in the OBO whose every target
            was obsolete and unresolvable (so they are absent from
            ``source_to_doids``).
    """

    xref_prefix: str
    source_to_doids: Dict[str, Set[str]]
    n_terms: int = 0
    n_obsolete_terms: int = 0
    n_xrefs: int = 0
    n_source_ids: int = 0
    n_non_numeric: int = 0
    n_one_to_many: int = 0
    n_obsolete_resolved: int = 0
    n_obsolete_dropped: int = 0
    source_ids_without_target: FrozenSet[str] = frozenset()

    def __len__(self) -> int:
        return len(self.source_to_doids)

    def targets(self, source_id: str) -> Set[str]:
        """Live DO terms cross-referencing ``source_id`` (empty if unmapped)."""
        return self.source_to_doids.get(source_id, set())


def _resolve_replacement(
    term: str, replaced_by: Dict[str, str], obsolete: Set[str]
) -> Optional[str]:
    """Follow ``replaced_by`` from an obsolete term to a live one.

    Chains are followed transitively (DO occasionally obsoletes a term that had
    itself replaced another) with a cycle guard. Returns ``None`` when the chain
    ends on an obsolete term, a missing term, or a cycle.
    """
    seen: Set[str] = {term}
    current = term
    while current in obsolete:
        nxt = replaced_by.get(current)
        if nxt is None or nxt in seen:
            return None
        seen.add(nxt)
        current = nxt
    return current


def build_doid_xref_map(path: Path, xref_prefix: str = "MIM") -> XrefMapping:
    """Build ``{external id → {DOID}}`` from the ``xref`` lines of ``doid.obo``.

    Args:
        path: the Disease Ontology OBO file.
        xref_prefix: which cross-reference namespace to index. ``"MIM"`` is
            OMIM (what UniProt ``DR MIM`` carries — the prefix in the OBO is
            ``MIM``, *not* ``OMIM``); ``"ORDO"`` is Orphanet.

    Returns:
        an :class:`XrefMapping` whose targets are live DO terms only, plus the
        counts needed to audit what the mapping dropped or duplicated.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Disease Ontology OBO not found: {path}")

    # Pass 1: collect every stanza's id, obsolete flag, replacement and xrefs.
    raw_xrefs: Dict[str, Set[str]] = defaultdict(set)  # source id → DO terms
    obsolete: Set[str] = set()
    replaced_by: Dict[str, str] = {}
    n_terms = 0
    n_xrefs = 0

    term_id: Optional[str] = None
    in_term = False

    with open(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("["):
                in_term = line == "[Term]"
                term_id = None
                if in_term:
                    n_terms += 1
            elif not in_term or not line:
                continue
            elif (match := _ID_RE.match(line)) is not None:
                term_id = match.group(1)
            elif term_id is None:
                continue
            elif line.startswith("is_obsolete:"):
                if line.split(":", 1)[1].strip().lower() == "true":
                    obsolete.add(term_id)
            elif (match := _REPLACED_BY_RE.match(line)) is not None:
                replaced_by[term_id] = match.group(1)
            elif (match := _XREF_RE.match(line)) is not None:
                if match.group(1) == xref_prefix:
                    raw_xrefs[match.group(2)].add(term_id)
                    n_xrefs += 1

    # Pass 2: resolve obsolete targets and record what each decision cost.
    source_to_doids: Dict[str, Set[str]] = {}
    n_obsolete_resolved = 0
    n_obsolete_dropped = 0
    without_target: Set[str] = set()

    for source_id, targets in raw_xrefs.items():
        live: Set[str] = set()
        for target in targets:
            if target not in obsolete:
                live.add(target)
                continue
            replacement = _resolve_replacement(target, replaced_by, obsolete)
            if replacement is None:
                n_obsolete_dropped += 1
            else:
                n_obsolete_resolved += 1
                live.add(replacement)
        if live:
            source_to_doids[source_id] = live
        else:
            without_target.add(source_id)

    n_non_numeric = sum(1 for source_id in raw_xrefs if not source_id.isdigit())
    n_one_to_many = sum(1 for targets in source_to_doids.values() if len(targets) > 1)

    mapping = XrefMapping(
        xref_prefix=xref_prefix,
        source_to_doids=source_to_doids,
        n_terms=n_terms,
        n_obsolete_terms=len(obsolete),
        n_xrefs=n_xrefs,
        n_source_ids=len(raw_xrefs),
        n_non_numeric=n_non_numeric,
        n_one_to_many=n_one_to_many,
        n_obsolete_resolved=n_obsolete_resolved,
        n_obsolete_dropped=n_obsolete_dropped,
        source_ids_without_target=frozenset(without_target),
    )
    logger.info(
        f"Parsed Disease Ontology {xref_prefix} cross-references from {path}: "
        f"{mapping.n_terms:,} terms ({mapping.n_obsolete_terms:,} obsolete), "
        f"{mapping.n_xrefs:,} {xref_prefix} xrefs over {mapping.n_source_ids:,} "
        f"distinct ids → {len(mapping):,} mappable"
    )
    logger.info(
        f"  one-to-many: {mapping.n_one_to_many:,} ids map to >1 DO term (kept, "
        f"expanded); non-numeric (phenotypic series): {mapping.n_non_numeric:,}; "
        f"obsolete targets: {mapping.n_obsolete_resolved:,} resolved via "
        f"replaced_by, {mapping.n_obsolete_dropped:,} dropped"
    )
    return mapping


def remap_protein_terms(
    protein_terms: Dict[str, Set[str]],
    mapping: XrefMapping,
    label: str = "OMIM→DOID",
) -> Tuple[Dict[str, Set[str]], RemapCoverage]:
    """Re-key a ``{protein: {source id}}`` map onto DO terms.

    One-to-many source ids expand to every target; unmapped ones are dropped
    (and counted); proteins left with no term at all disappear from the map, as
    they must — the Fisher engine treats an absent protein as having no
    annotation, which is exactly right for a protein whose only diseases are
    outside DO.

    A thin DOID-flavoured wrapper over the generic :func:`src.remap.remap_values`
    (which the gene-keyed layers also use, on the protein axis).

    Returns:
        ``(remapped map, coverage)``. The coverage report is returned rather than
        only logged so callers and tests can assert on it.
    """
    return remap_values(
        protein_terms,
        mapping,
        label,
        key_label="protein",
        value_label="term",
        target_label="Disease Ontology term",
    )


class DiseaseOntologyAnnotationSource(AnnotationSource):
    """Domain annotations keyed by Disease Ontology term.

    Reads a UniProt ``DR`` disease vocabulary (OMIM phenotypes by default,
    Orphanet with ``xref_prefix="ORDO"``) and re-keys it to DOID terms *before*
    the annotations reach the statistics, so sparse per-locus disease entries
    pool into the DO class that subsumes them.
    """

    def __init__(
        self,
        dat_path: Path,
        doid_obo_path: Path,
        database: str = "MIM",
        id_type: Optional[str] = "phenotype",
        xref_prefix: str = "MIM",
        spec: OntologySpec = DOID_SPEC,
    ) -> None:
        self.dat_path = Path(dat_path)
        self.doid_obo_path = Path(doid_obo_path)
        self.database = database
        self.id_type = id_type
        self.xref_prefix = xref_prefix
        self.spec = spec
        #: Populated by :meth:`parse`, for callers that want the audit numbers.
        self.coverage: Optional[RemapCoverage] = None

    def parse(self) -> Dict[str, Set[str]]:
        raw = parse_uniprot_cross_refs(
            self.dat_path, self.database, id_type=self.id_type
        )
        mapping = build_doid_xref_map(self.doid_obo_path, self.xref_prefix)
        remapped, self.coverage = remap_protein_terms(
            raw, mapping, label=f"{self.xref_prefix}→DOID"
        )
        return remapped
