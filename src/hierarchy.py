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
"""

from __future__ import annotations

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
