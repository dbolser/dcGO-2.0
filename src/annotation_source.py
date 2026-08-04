"""Annotation-source abstraction: the seam for associating domains with any ontology.

The dcGO statistical engine (``sparse_fisher`` → ``vectorized_fisher``) is
completely agnostic to what an "annotation term" means. It consumes two plain
dictionaries — ``{protein_id: {domain}}`` and ``{protein_id: {term}}`` — and
never inspects the term namespace. All the coupling to Gene Ontology lives at
the *input boundary*: how a protein→term map is built.

An :class:`AnnotationSource` is that boundary. Each concrete source parses one
kind of annotation file and returns ``{protein_id: {term_id}}`` keyed by the
**same protein id space as the domain annotations** (UniProt accessions, as in
``protein2ipr``). Adding a new ontology (Disease Ontology, HPO, EC, Reactome …)
means writing a new ``AnnotationSource`` subclass plus, for hierarchical
ontologies, pointing :class:`OntologySpec` at its DAG — the Fisher/FDR core and
the True Path propagation are reused unchanged.

The reference implementation, :class:`GAFAnnotationSource`, wraps the existing
GOA parser so the current human GO path routes through this seam.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Set


@dataclass(frozen=True)
class OntologySpec:
    """Metadata describing an ontology layer the pipeline can associate with domains.

    Attributes:
        ontology_id: Short stable identifier, e.g. ``"GO"``, ``"DO"``, ``"HPO"``, ``"EC"``.
        name: Human-readable name, e.g. ``"Gene Ontology"``.
        term_prefix: Prefix every term id carries (e.g. ``"GO:"``), used to
            validate/recognise this ontology's terms. ``None`` = no prefix check
            (e.g. EC numbers like ``1.1.1.1``).
        obo_path: Path to the ontology DAG in OBO format, for True Path Rule
            propagation. ``None`` for ontologies without an OBO graph (a
            non-hierarchical layer, or one needing a bespoke hierarchy loader
            such as EC/Reactome).
    """

    ontology_id: str
    name: str
    term_prefix: Optional[str] = None
    obo_path: Optional[Path] = None


# The default spec for Gene Ontology — the ontology the pipeline ships with.
GO_SPEC = OntologySpec(
    ontology_id="GO",
    name="Gene Ontology",
    term_prefix="GO:",
)


class AnnotationSource(ABC):
    """A source of ``{protein_id: {ontology_term}}`` annotations for one ontology.

    Concrete subclasses own the parsing and — critically — any identifier
    mapping needed to key the result by the same protein id space as the domain
    annotations (UniProt accessions). Ontologies whose annotations arrive keyed
    by gene (HGNC/OMIM/Ensembl/MGI) must resolve to UniProt inside ``parse`` so
    the engine can join domains and terms on a shared key.
    """

    #: Metadata for the ontology this source provides. Subclasses must set this.
    spec: OntologySpec

    @abstractmethod
    def parse(self) -> Dict[str, Set[str]]:
        """Return a mapping of protein id → set of ontology term ids."""
        raise NotImplementedError


class GAFAnnotationSource(AnnotationSource):
    """Gene Ontology annotations from a GAF 2.2 file (GOA).

    The reference :class:`AnnotationSource`: it wraps :func:`src.goa_parser.parse_goa`,
    which is species-agnostic (``goa_human.gaf.gz``, ``goa_mouse.gaf.gz``, …).
    """

    def __init__(
        self,
        gaf_path: Path,
        evidence_filter: str = "manual",
        aspects: Optional[Set[str]] = None,
        spec: OntologySpec = GO_SPEC,
    ) -> None:
        self.gaf_path = Path(gaf_path)
        self.evidence_filter = evidence_filter
        self.aspects = aspects or {"P", "F", "C"}
        self.spec = spec

    def parse(self) -> Dict[str, Set[str]]:
        # Imported lazily so the abstraction module stays dependency-light.
        from src.goa_parser import parse_goa

        return parse_goa(
            self.gaf_path,
            evidence_filter=self.evidence_filter,
            aspects=self.aspects,
        )


def restrict_to_universe(
    protein_terms: Dict[str, Set[str]], universe: Set[str]
) -> Dict[str, Set[str]]:
    """Drop annotations for proteins outside the analysis universe.

    The Fisher engine's protein universe is the *intersection* of the domain map
    and the annotation map — a protein with no domain assignment is missing
    data, not evidence that a domain is absent. ``build_sparse_matrices`` keys
    its rows off the **union** of the two maps, so the caller has to narrow the
    annotation map itself.

    This matters far more for the UniProt-native sources than it does for GOA.
    ``protein2ipr`` is extracted per species, but ``uniprot_sprot.dat`` is the
    whole of Swiss-Prot, so a ``--species human --ontology reactome`` run parsed
    39,418 Reactome-annotated proteins from every organism; without this
    restriction each of the ~28k non-human ones entered every contingency table
    as a domain-negative observation and shifted the term backgrounds.
    """
    return {
        protein: terms
        for protein, terms in protein_terms.items()
        if protein in universe
    }
