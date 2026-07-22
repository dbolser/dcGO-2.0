"""Enzyme Commission (EC) annotations as an :class:`AnnotationSource`.

EC is the simplest non-GO ontology to plug into the dcGO engine because it
sidesteps the identifier-mapping problem: the Expasy **ENZYME** database
(``enzyme.dat``) already lists, for each EC number, the UniProt accessions
annotated with it (``DR`` lines) — the *same* id space as ``protein2ipr``. So
``parse`` yields ``{uniprot_accession: {ec_number}}`` that joins directly to the
domain annotations, with no HGNC/Ensembl reconciliation.

The EC hierarchy is implicit in the dotted numbering
(``1.1.1.1`` ⊂ ``1.1.1.-`` ⊂ ``1.1.-.-`` ⊂ ``1.-.-.-``), so ancestor lookup is
pure string manipulation (:func:`ec_ancestors`) — no OBO graph required. That
makes a True-Path-style propagation for EC a small, self-contained follow-up
rather than an obonet dependency.

enzyme.dat record shape (Swiss-Prot flat-file style, ``//``-delimited)::

    ID   1.1.1.1
    DE   Alcohol dehydrogenase.
    DR   P07327, ADH1A_HUMAN;  P28469, ADH1A_MACMU;  Q5RBP7, ADH1A_PONAB;
    //

"Transferred entry" and "Deleted entry" records carry no ``DR`` lines and are
ignored.
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, Set

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec

if TYPE_CHECKING:
    from src.ontology_processor import Annotation

# EC numbers carry no id prefix (e.g. "1.1.1.1"), so term_prefix stays None.
EC_SPEC = OntologySpec(ontology_id="EC", name="Enzyme Commission")

_EC_LEVELS = 4


def ec_ancestors(ec_number: str) -> List[str]:
    """Return the ancestor EC numbers of ``ec_number``, most specific first.

    ``"1.1.1.1"`` → ``["1.1.1.-", "1.1.-.-", "1.-.-.-"]``. Partial numbers are
    handled: ``"1.1.1.-"`` → ``["1.1.-.-", "1.-.-.-"]``. Anything that is not a
    4-level dotted number (e.g. a malformed id) yields ``[]``.
    """
    parts = ec_number.split(".")
    if len(parts) != _EC_LEVELS:
        return []

    ancestors: List[str] = []
    for i in range(_EC_LEVELS - 1, 0, -1):
        if parts[i] == "-":
            continue
        parent = parts.copy()
        for j in range(i, _EC_LEVELS):
            parent[j] = "-"
        ancestors.append(".".join(parent))
    return ancestors


def propagate_ec_annotations(direct_associations: Iterable[Any]) -> "List[Annotation]":
    """Propagate significant domain→EC associations up the EC hierarchy.

    The EC counterpart of the GO True Path Rule
    (``ontology_processor.propagate_annotations``): instead of walking an OBO
    DAG, it uses the implicit EC hierarchy via :func:`ec_ancestors`. Each direct
    association is emitted as-is, plus one ``"propagated"`` annotation per
    ancestor EC number. ``(domain, ec)`` pairs are de-duplicated, keeping the
    most significant source association (lowest ``q_value``), so a shared
    ancestor is attributed to the strongest evidence.

    Args:
        direct_associations: iterable of objects exposing ``domain``,
            ``go_term`` (the EC number), ``q_value`` and ``hyper_score`` — e.g.
            the ``AssociationResult`` records built in ``run_dcgo_human.py``.

    Returns:
        list of ``ontology_processor.Annotation`` with ``annotation_type`` of
        ``"direct"`` or ``"propagated"``.
    """
    # Lazy import keeps this module dependency-light (ontology_processor pulls in
    # obonet/networkx/pandas) and avoids an import cycle.
    from src.ontology_processor import Annotation

    # Most significant first, so shared ancestors record the best source term.
    ordered = sorted(direct_associations, key=lambda a: a.q_value)

    annotations = []
    seen: Set[tuple] = set()  # (domain, ec) pairs already emitted
    for assoc in ordered:
        direct_key = (assoc.domain, assoc.go_term)
        if direct_key not in seen:
            annotations.append(
                Annotation(
                    domain=assoc.domain,
                    go_term=assoc.go_term,
                    q_value=assoc.q_value,
                    association_score=assoc.hyper_score,
                    annotation_type="direct",
                    direct_source_term=assoc.go_term,
                )
            )
            seen.add(direct_key)

        for ancestor in ec_ancestors(assoc.go_term):
            ancestor_key = (assoc.domain, ancestor)
            if ancestor_key not in seen:
                annotations.append(
                    Annotation(
                        domain=assoc.domain,
                        go_term=ancestor,
                        q_value=assoc.q_value,
                        association_score=assoc.hyper_score,
                        annotation_type="propagated",
                        direct_source_term=assoc.go_term,
                    )
                )
                seen.add(ancestor_key)

    return annotations


def parse_enzyme_dat(path: Path) -> Dict[str, Set[str]]:
    """Parse Expasy ``enzyme.dat`` into ``{uniprot_accession: {ec_number}}``.

    Args:
        path: Path to ``enzyme.dat`` (optionally gzipped).

    Returns:
        Mapping of UniProt accession → set of EC numbers it is annotated with.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"ENZYME data file not found: {path}")

    logger.info(f"Parsing ENZYME database: {path}")
    open_func = gzip.open if path.suffix == ".gz" else open

    protein_ec: Dict[str, Set[str]] = defaultdict(set)
    current_ec: str | None = None
    n_active = 0
    n_skipped = 0

    with open_func(path, "rt") as f:
        for line in f:
            tag = line[:2]

            if tag == "ID":
                current_ec = line[5:].strip() or None
            elif tag == "DE" and current_ec is not None:
                desc = line[5:].strip()
                # Transferred/deleted entries are placeholders, not real enzymes.
                if desc.startswith("Transferred entry") or desc.startswith(
                    "Deleted entry"
                ):
                    current_ec = None
                    n_skipped += 1
            elif tag == "DR" and current_ec is not None:
                # "P07327, ADH1A_HUMAN;  P28469, ADH1A_MACMU;" → first field of
                # each ';'-separated pair is the UniProt accession.
                for chunk in line[5:].split(";"):
                    chunk = chunk.strip()
                    if not chunk:
                        continue
                    accession = chunk.split(",")[0].strip()
                    if accession:
                        protein_ec[accession].add(current_ec)
            elif line.startswith("//"):
                if current_ec is not None:
                    n_active += 1
                current_ec = None

    all_ec = set().union(*protein_ec.values()) if protein_ec else set()
    logger.info("ENZYME parsing complete:")
    logger.info(f"  Active EC entries with proteins: {n_active:,}")
    logger.info(f"  Transferred/deleted entries skipped: {n_skipped:,}")
    logger.info(f"  Unique enzymes (proteins): {len(protein_ec):,}")
    logger.info(f"  Unique EC numbers: {len(all_ec):,}")

    return dict(protein_ec)


class ECAnnotationSource(AnnotationSource):
    """Enzyme Commission annotations from the Expasy ENZYME database.

    A concrete :class:`AnnotationSource`: the dcGO engine associates protein
    domains with EC numbers exactly as it does with GO terms, since both arrive
    as ``{protein_id: {term_id}}`` keyed by UniProt accession.
    """

    def __init__(self, enzyme_dat_path: Path, spec: OntologySpec = EC_SPEC) -> None:
        self.enzyme_dat_path = Path(enzyme_dat_path)
        self.spec = spec

    def parse(self) -> Dict[str, Set[str]]:
        return parse_enzyme_dat(self.enzyme_dat_path)
