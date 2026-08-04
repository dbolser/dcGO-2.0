"""Generic term-hierarchy propagation (the True Path Rule, ontology-agnostic).

The dcGO engine produces *direct* domain→term associations. The True Path Rule
then propagates each up its ontology's hierarchy so a domain credited with a
specific term is also credited with the more general ancestors. GO does this via
an OBO DAG (``ontology_processor``); EC via its dotted numbering
(``ec_ancestors``); Reactome via its pathway relations; UniProt keywords via
their keyword DAG.

All of those differ only in *how you find a term's ancestors*. This module
captures the shared machinery:

* :func:`closure_ancestors` turns a ``{child: {parents}}`` map (Reactome
  relations, keyword hierarchy, …) into a memoised transitive-ancestors function.
* :func:`propagate_via_ancestors` takes the direct associations and any
  ancestors function and emits direct + propagated annotations, de-duplicating
  ``(domain, term)`` pairs and attributing each to its most significant source.

Two families of hierarchy *loaders* live here too, because they are ontology
agnostic:

* **Implicit hierarchies encoded in the id itself** — :func:`dotted_ancestors`
  (EC ``1.1.1.1``, TCDB ``8.A.98.1.10``) and :func:`alpha_prefix_ancestors`
  (MEROPS ``S01.151`` → ``S01`` → ``S``, CAZy ``GT32`` → ``GT``).
* **OBO graphs** — :func:`parse_obo_child_parents` reads ``is_a`` (and,
  optionally, ``part_of``) edges out of any OBO file into a child→parents map,
  which is all ChEBI/DO/HPO-style ontologies need to reach
  :func:`closure_ancestors`. (GO keeps its richer obonet-based
  ``ontology_processor`` path, which also does optimal-level filtering.)
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, Iterable, List, Set

if TYPE_CHECKING:
    from src.ontology_processor import Annotation


def closure_ancestors(
    child_to_parents: Dict[str, Set[str]],
) -> Callable[[str], Set[str]]:
    """Build a memoised transitive-ancestors function from a child→parents map.

    ``ancestors(term)`` returns every term reachable by walking ``parents``
    edges upward (excluding ``term`` itself). Safe on DAGs; a cycle guard keeps
    accidental cycles from recursing forever (the offending back-edge is simply
    not expanded).
    """
    cache: Dict[str, Set[str]] = {}

    def ancestors(term: str) -> Set[str]:
        if term in cache:
            return cache[term]
        cache[term] = set()  # cycle guard: a re-entrant visit sees an empty set
        result: Set[str] = set()
        for parent in child_to_parents.get(term, ()):
            result.add(parent)
            result |= ancestors(parent)
        cache[term] = result
        return result

    return ancestors


def propagate_via_ancestors(
    direct_associations: Iterable[Any],
    ancestors_fn: Callable[[str], Iterable[str]],
) -> "List[Annotation]":
    """Propagate significant domain→term associations up a term hierarchy.

    Each direct association is emitted as-is, plus one ``"propagated"``
    annotation per ancestor term returned by ``ancestors_fn``. ``(domain, term)``
    pairs are de-duplicated, keeping the most significant source association
    (lowest ``q_value``), so a shared ancestor is attributed to the strongest
    evidence.

    Args:
        direct_associations: iterable of objects exposing ``domain``,
            ``go_term`` (the term id, whatever the ontology), ``q_value`` and
            ``hyper_score`` — e.g. the ``AssociationResult`` records built in
            ``run_dcgo_human.py``.
        ancestors_fn: maps a term id to its ancestor term ids.

    Returns:
        list of ``ontology_processor.Annotation`` (``annotation_type`` of
        ``"direct"`` or ``"propagated"``).
    """
    # Lazy import keeps this module dependency-light (ontology_processor pulls in
    # obonet/networkx/pandas) and avoids an import cycle.
    from src.ontology_processor import Annotation

    # Most significant first, so shared ancestors record the best source term.
    ordered = sorted(direct_associations, key=lambda a: a.q_value)

    annotations: List[Annotation] = []
    seen: Set[tuple] = set()  # (domain, term) pairs already emitted
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

        for ancestor in ancestors_fn(assoc.go_term):
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


def dotted_ancestors(term: str, separator: str = ".") -> List[str]:
    """Ancestors of a dotted classification id, most specific first.

    For classifications whose hierarchy *is* the id — TCDB
    (``8.A.98.1.10`` ⊂ ``8.A.98.1`` ⊂ … ⊂ ``8``), Pathway/Brite-style numbering,
    Rhea-free EC variants — each ancestor is the id truncated by one level::

        >>> dotted_ancestors("8.A.98.1.10")
        ['8.A.98.1', '8.A.98', '8.A', '8']

    A single-level id has no ancestors. (EC keeps its own
    :func:`src.ec_annotation_source.ec_ancestors`, which pads with ``-``
    placeholders instead of truncating, because that is how EC writes its
    partial numbers.)
    """
    parts = [p for p in term.split(separator) if p]
    if len(parts) < 2:
        return []
    return [separator.join(parts[:i]) for i in range(len(parts) - 1, 0, -1)]


def alpha_prefix_ancestors(term: str) -> List[str]:
    """Ancestors of an id shaped ``<letters><digits>[.<subid>]``, most specific first.

    Covers the family/class classifications whose hierarchy is spelled out in
    the accession: MEROPS peptidases (``S01.151`` → family ``S01`` → catalytic
    type ``S``) and CAZy families (``GT32`` → class ``GT``)::

        >>> alpha_prefix_ancestors("S01.151")
        ['S01', 'S']
        >>> alpha_prefix_ancestors("GT32")
        ['GT']

    Ids that do not start with letters followed by digits yield ``[]``.
    """
    match = re.match(r"^([A-Za-z]+)(\d+)", term)
    if not match:
        return []
    letters, digits = match.group(1), match.group(2)
    family = f"{letters}{digits}"
    ancestors: List[str] = []
    if term != family:  # e.g. "S01.151" → its family "S01"
        ancestors.append(family)
    ancestors.append(letters)
    return ancestors


def parse_obo_child_parents(
    path: Path,
    relations: Iterable[str] = ("part_of",),
    include_obsolete: bool = False,
) -> Dict[str, Set[str]]:
    """Read an OBO file into ``{child_id: {parent_ids}}``.

    A deliberately small OBO reader: it walks ``[Term]`` stanzas and keeps
    ``is_a`` edges plus any ``relationship: <rel> <id>`` whose ``<rel>`` is in
    ``relations`` (``part_of`` by default — the other relation the True Path
    Rule traverses). That is everything :func:`closure_ancestors` needs, and it
    avoids pulling obonet/networkx in for ontologies as large as ChEBI.

    Args:
        path: OBO file (optionally gzipped).
        relations: non-``is_a`` relationship types to treat as parent edges.
        include_obsolete: keep stanzas marked ``is_obsolete: true``.

    Returns:
        child → set of direct parents.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"OBO file not found: {path}")

    wanted = set(relations)
    open_func = gzip.open if path.suffix == ".gz" else open

    child_to_parents: Dict[str, Set[str]] = defaultdict(set)
    term_id: str | None = None
    parents: Set[str] = set()
    obsolete = False
    in_term = False

    def flush() -> None:
        if term_id and (include_obsolete or not obsolete) and parents:
            child_to_parents[term_id] |= parents

    with open_func(path, "rt") as f:
        for raw in f:
            line = raw.strip()
            if line.startswith("["):
                flush()
                in_term = line == "[Term]"
                term_id, parents, obsolete = None, set(), False
            elif not in_term or not line:
                continue
            elif line.startswith("id:"):
                term_id = line[3:].strip()
            elif line.startswith("is_a:"):
                # "is_a: CHEBI:24431 ! chemical entity"
                parents.add(line[5:].split("!")[0].strip())
            elif line.startswith("is_obsolete:"):
                obsolete = line.split(":", 1)[1].strip().lower() == "true"
            elif line.startswith("relationship:"):
                fields = line[13:].split("!")[0].split()
                if len(fields) >= 2 and fields[0] in wanted:
                    parents.add(fields[1].strip())
        flush()

    return dict(child_to_parents)
