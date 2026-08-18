"""CIViC clinical evidence re-keyed onto NCIt and OncoTree cancer vocabularies.

CIViC (CC0) curates clinical evidence for cancer variants: each evidence row
names a *molecular profile* (a variant of one or more genes) and a disease,
DOID-keyed. Collapsing the variant axis gives gene → cancer-type annotations:
"variants of this gene have clinical significance in this cancer". Neither
NCIt nor OncoTree publishes a direct gene → term file, so this chain —
CIViC gene → DOID → {NCIt, OncoTree} via cross-reference tables — is the best
honest route, and every stage of it is counted.

Row policy (:class:`CIViCFilterCounts`):

* Only the evidence summaries file is read (4,886 rows at the 2026-08 nightly;
  the ~150 assertions are curated *over the same evidence items*, so adding
  them would double-count evidence without adding genes).
* ``evidence_direction`` must be ``Supports`` — "Does Not Support" rows (518)
  are negative evidence and ``N/A`` (17) is no direction at all; both dropped
  and counted. The nightly ships only ``accepted`` rows, but the status is
  checked and counted anyway.
* Rows with no DOID (230) are dropped and counted: the disease cannot enter
  either target vocabulary without one.

**Gene extraction** is the one heuristic step: CIViC names profiles
``<GENE> <variant>`` (``JAK2 V617F``), joins complex profiles with ``AND`` /
``OR`` (``BCR::ABL1 Fusion AND ABL1 T315I``) and writes fusions
``A::B``. Each component's first whitespace token is taken as its gene, and
fusion tokens credit both partners. These names are machine-generated from
CIViC's gene records, so the first token is reliable; anything that is not a
real symbol simply fails the audited symbol → UniProt re-key later rather
than passing silently.

The disease axis then re-keys DOID → target, and the audit matters because
the two targets differ sharply (at DO 2026-07-31 / OncoTree 2025-10-03):

* **NCIt** — ``doid.obo`` carries ``xref: NCI:C…`` on 4,737 terms; 247 of the
  299 CIViC DOIDs map, reaching 455 of 497 genes.
* **OncoTree** — no DO xref exists, so the chain goes through the OncoTree
  API's own ``externalReferences`` (NCI and UMLS ids per node), matched
  against DO's ``NCI`` and ``UMLS_CUI`` xrefs. Only 124 of 299 DOIDs map
  (310 genes): OncoTree is a clinical *tumour-type* tree and simply has no
  node for many DO classes CIViC uses (gene-level "cancer" umbrella terms,
  benign conditions). That thinness is inherent to the target, not a defect
  of the chain, and the per-run coverage log states it.

Both layers finish with the standard audited gene-symbol → UniProt re-key
(Swiss-Prot ``DR HGNC`` symbols). Hierarchies: ``ncit.obo`` (``is_a``; the
shared OBO reader streams it, so its 248 MB never sits in memory — only the
child→parents map does) and the OncoTree JSON's ``parent`` field (every node
chains up to OncoTree's ``TISSUE`` root).
"""

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec
from src.disease_ontology import build_doid_xref_map, build_doid_xref_maps
from src.gene_mapping import parse_gene_accession_index, remap_gene_annotations
from src.remap import RemapCoverage, remap_values

#: NCI Thesaurus terms, e.g. ``NCIT:C3058``.
NCIT_SPEC = OntologySpec(ontology_id="NCIT", name="NCI Thesaurus", term_prefix="NCIT:")

#: OncoTree tumour-type codes, e.g. ``GBM`` (flat codes; the tree is the
#: JSON's parent field).
ONCOTREE_SPEC = OntologySpec(ontology_id="OncoTree", name="OncoTree tumour type")

_PROFILE_SPLIT_RE = re.compile(r"\s+(?:AND|OR)\s+")


@dataclass(frozen=True)
class CIViCFilterCounts:
    """What the evidence-row policy did.

    Attributes:
        n_rows: evidence rows read.
        n_not_supports: rows whose direction was not ``Supports``.
        n_not_accepted: rows whose status was not ``accepted`` (0 in the
            nightly, counted defensively).
        n_no_doid: rows with an empty DOID column.
        n_no_gene: rows whose molecular-profile name yielded no gene token —
            its own counter so the audit names the stage that dropped the row.
        n_kept: rows contributing at least one (gene, DOID) pair.
    """

    n_rows: int = 0
    n_not_supports: int = 0
    n_not_accepted: int = 0
    n_no_doid: int = 0
    n_no_gene: int = 0
    n_kept: int = 0


def genes_of_molecular_profile(profile: str) -> Set[str]:
    """The gene symbols a CIViC molecular-profile name involves.

    ``"JAK2 V617F"`` → ``{"JAK2"}``; ``"BCR::ABL1 Fusion AND ABL1 T315I"`` →
    ``{"BCR", "ABL1"}``. See the module docstring for why the first token of
    each AND/OR component is the gene.
    """
    genes: Set[str] = set()
    for component in _PROFILE_SPLIT_RE.split(profile.strip()):
        tokens = component.split()
        if not tokens:
            continue
        for gene in tokens[0].split("::"):
            if gene:
                genes.add(gene)
    return genes


def parse_civic_evidence(
    path: Path,
) -> Tuple[Dict[str, Set[str]], CIViCFilterCounts]:
    """Parse CIViC evidence summaries into ``{gene symbol: {DOID term}}``."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"CIViC evidence summaries file not found: {path}")

    logger.info(f"Parsing CIViC clinical evidence from {path}")
    gene_terms: Dict[str, Set[str]] = defaultdict(set)
    n_rows = n_not_supports = n_not_accepted = n_no_doid = n_no_gene = n_kept = 0
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"molecular_profile", "doid", "evidence_direction"} - set(
            reader.fieldnames or ()
        )
        if missing:
            raise ValueError(
                f"{path} lacks expected column(s) {sorted(missing)}; "
                f"header was {reader.fieldnames}"
            )
        for row in reader:
            n_rows += 1
            if (row.get("evidence_status") or "accepted").strip() != "accepted":
                n_not_accepted += 1
                continue
            if (row["evidence_direction"] or "").strip() != "Supports":
                n_not_supports += 1
                continue
            doid = (row["doid"] or "").strip()
            if not doid:
                n_no_doid += 1
                continue
            genes = genes_of_molecular_profile(row["molecular_profile"] or "")
            if not genes:
                n_no_gene += 1
                continue
            n_kept += 1
            for gene in genes:
                gene_terms[gene].add(f"DOID:{doid}")

    counts = CIViCFilterCounts(
        n_rows=n_rows,
        n_not_supports=n_not_supports,
        n_not_accepted=n_not_accepted,
        n_no_doid=n_no_doid,
        n_no_gene=n_no_gene,
        n_kept=n_kept,
    )
    logger.info(
        f"  Rows: {counts.n_rows:,}; not 'Supports': {counts.n_not_supports:,}; "
        f"not accepted: {counts.n_not_accepted:,}; no DOID: "
        f"{counts.n_no_doid:,}; no gene token: {counts.n_no_gene:,}; "
        f"kept: {counts.n_kept:,}"
    )
    logger.info(
        f"  Genes: {len(gene_terms):,}; diseases: "
        f"{len({term for terms in gene_terms.values() for term in terms}):,}"
    )
    return dict(gene_terms), counts


@dataclass(frozen=True)
class TermTargetMap:
    """A plain ``source term → target terms`` table for :func:`remap_values`."""

    source_to_targets: Dict[str, Set[str]]

    def __len__(self) -> int:
        return len(self.source_to_targets)

    def targets(self, source_id: str) -> Set[str]:
        return self.source_to_targets.get(source_id, set())


def build_doid_to_ncit_map(doid_obo_path: Path) -> TermTargetMap:
    """``{DOID term → {NCIT term}}`` from doid.obo's ``NCI`` cross-references.

    Reuses the audited DO xref reader (which resolves obsolete DO stanzas
    through ``replaced_by``) and inverts its NCI-id → DOID direction; the NCI
    ids gain the ``NCIT:`` prefix ncit.obo uses.
    """
    mapping = build_doid_xref_map(doid_obo_path, "NCI")
    doid_to_ncit: Dict[str, Set[str]] = defaultdict(set)
    for nci_id, doids in mapping.source_to_doids.items():
        for doid in doids:
            doid_to_ncit[doid].add(f"NCIT:{nci_id}")
    logger.info(
        f"  Inverted: {len(doid_to_ncit):,} DOID terms carry an NCI cross-reference"
    )
    return TermTargetMap(dict(doid_to_ncit))


@dataclass(frozen=True)
class OncoTreeVocabulary:
    """The OncoTree tumour-type tree and its external references.

    Attributes:
        child_to_parents: ``code → {parent code}`` (roots chain to ``TISSUE``).
        nci_to_codes / umls_to_codes: reverse indexes over each node's
            ``externalReferences``.
        n_nodes: tumour types read.
    """

    child_to_parents: Dict[str, Set[str]]
    nci_to_codes: Dict[str, Set[str]]
    umls_to_codes: Dict[str, Set[str]]
    n_nodes: int = 0


@lru_cache(maxsize=8)
def _load_oncotree(path_str: str) -> OncoTreeVocabulary:
    """Parse and validate one OncoTree dump; memoised so the ``oncotree``
    entry's annotation chain and hierarchy share a single parse per file."""
    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(f"OncoTree tumour-types JSON not found: {path}")
    with open(path, "rt", encoding="utf-8") as handle:
        nodes = json.load(handle)
    # The API can serve an error/envelope object instead of the node list;
    # accepting it would AttributeError on a string node or, worse, silently
    # build an empty vocabulary. Validate the shape and name the file.
    if not isinstance(nodes, list) or not nodes:
        raise ValueError(
            f"{path} does not look like the OncoTree tumorTypes dump: expected "
            f"a non-empty JSON list of node objects, got "
            f"{type(nodes).__name__}"
        )
    bad = next((node for node in nodes if not isinstance(node, dict)), None)
    if bad is not None:
        raise ValueError(
            f"{path} does not look like the OncoTree tumorTypes dump: list "
            f"entries should be node objects, found {type(bad).__name__}"
        )
    if not any("code" in node for node in nodes):
        raise ValueError(
            f"{path} does not look like the OncoTree tumorTypes dump: no "
            "entry carries a 'code' field"
        )
    child_to_parents: Dict[str, Set[str]] = {}
    nci_to_codes: Dict[str, Set[str]] = defaultdict(set)
    umls_to_codes: Dict[str, Set[str]] = defaultdict(set)
    for node in nodes:
        code = node.get("code")
        if not code:
            continue
        parent = node.get("parent")
        if parent:
            child_to_parents[code] = {parent}
        references = node.get("externalReferences") or {}
        for nci_id in references.get("NCI") or ():
            nci_to_codes[nci_id].add(code)
        for umls_id in references.get("UMLS") or ():
            umls_to_codes[umls_id].add(code)
    vocabulary = OncoTreeVocabulary(
        child_to_parents=child_to_parents,
        nci_to_codes=dict(nci_to_codes),
        umls_to_codes=dict(umls_to_codes),
        n_nodes=len(nodes),
    )
    logger.info(
        f"Parsed OncoTree from {path}: {vocabulary.n_nodes:,} tumour types, "
        f"{len(vocabulary.nci_to_codes):,} with NCI ids, "
        f"{len(vocabulary.umls_to_codes):,} with UMLS ids"
    )
    return vocabulary


def parse_oncotree(path: Path) -> OncoTreeVocabulary:
    """Parse the OncoTree API's tumour-types JSON (validated, memoised)."""
    return _load_oncotree(str(Path(path)))


def parse_oncotree_child_parents(path: Path) -> Dict[str, Set[str]]:
    """The OncoTree hierarchy alone, for the registry's ancestors/parents."""
    return parse_oncotree(path).child_to_parents


def build_doid_to_oncotree_map(
    doid_obo_path: Path, oncotree_path: Path
) -> TermTargetMap:
    """``{DOID term → {OncoTree code}}`` via NCI and UMLS cross-references.

    DO carries no OncoTree xref, so a DOID reaches a code when the two agree
    on an NCI Thesaurus id or a UMLS CUI (union of both routes; both xref
    namespaces come from one pass over doid.obo).
    """
    vocabulary = parse_oncotree(oncotree_path)
    mappings = build_doid_xref_maps(doid_obo_path, ("NCI", "UMLS_CUI"))
    doid_to_codes: Dict[str, Set[str]] = defaultdict(set)
    for prefix, reverse_index in (
        ("NCI", vocabulary.nci_to_codes),
        ("UMLS_CUI", vocabulary.umls_to_codes),
    ):
        for external_id, doids in mappings[prefix].source_to_doids.items():
            codes = reverse_index.get(external_id)
            if codes:
                for doid in doids:
                    doid_to_codes[doid] |= codes
    logger.info(
        f"  Chained: {len(doid_to_codes):,} DOID terms reach an OncoTree code "
        f"(via shared NCI/UMLS ids)"
    )
    return TermTargetMap(dict(doid_to_codes))


class _CIViCChainedAnnotationSource(AnnotationSource):
    """Shared chain: CIViC genes → DOID → target vocabulary → UniProt."""

    def __init__(self, civic_path: Path, dat_path: Path) -> None:
        self.civic_path = Path(civic_path)
        self.dat_path = Path(dat_path)
        #: Populated by :meth:`parse`: the evidence-row audit.
        self.filter_counts: Optional[CIViCFilterCounts] = None
        #: Populated by :meth:`parse`: the DOID → target stage (its *keys*
        #: are genes, its *values* the DOID terms being remapped).
        self.disease_coverage: Optional[RemapCoverage] = None
        #: Populated by :meth:`parse`: the gene → UniProt stage (axis-swapped
        #: as in remap_gene_annotations).
        self.coverage: Optional[RemapCoverage] = None

    def _build_disease_mapping(self) -> TermTargetMap:
        raise NotImplementedError

    def parse(self) -> Dict[str, Set[str]]:
        gene_doids, self.filter_counts = parse_civic_evidence(self.civic_path)
        mapping = self._build_disease_mapping()
        gene_terms, self.disease_coverage = remap_values(
            gene_doids,
            mapping,
            label=f"DOID→{self.spec.ontology_id} (CIViC)",
            key_label="gene",
            value_label="DOID",
            target_label=f"{self.spec.name} term",
        )
        index = parse_gene_accession_index(self.dat_path)
        remapped, self.coverage = remap_gene_annotations(
            gene_terms,
            index.symbol,
            label=f"gene symbol→UniProt (CIViC/{self.spec.ontology_id})",
        )
        return remapped


class CIViCNCItAnnotationSource(_CIViCChainedAnnotationSource):
    """Domain annotations keyed by NCI Thesaurus cancer term (via CIViC)."""

    spec = NCIT_SPEC

    def __init__(self, civic_path: Path, doid_obo_path: Path, dat_path: Path) -> None:
        super().__init__(civic_path, dat_path)
        self.doid_obo_path = Path(doid_obo_path)

    def _build_disease_mapping(self) -> TermTargetMap:
        return build_doid_to_ncit_map(self.doid_obo_path)


class CIViCOncoTreeAnnotationSource(_CIViCChainedAnnotationSource):
    """Domain annotations keyed by OncoTree tumour-type code (via CIViC)."""

    spec = ONCOTREE_SPEC

    def __init__(
        self,
        civic_path: Path,
        doid_obo_path: Path,
        oncotree_path: Path,
        dat_path: Path,
    ) -> None:
        super().__init__(civic_path, dat_path)
        self.doid_obo_path = Path(doid_obo_path)
        self.oncotree_path = Path(oncotree_path)

    def _build_disease_mapping(self) -> TermTargetMap:
        return build_doid_to_oncotree_map(self.doid_obo_path, self.oncotree_path)
