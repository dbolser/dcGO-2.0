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
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Sequence,
    Set,
)

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


def parents_from_map(
    child_to_parents: Dict[str, Set[str]],
) -> Callable[[str], Set[str]]:
    """Direct-parents function over a child→parents map.

    The map *is* the direct-parent relation, so this is a lookup — but relative
    inference needs it as a function, alongside the transitive
    :func:`closure_ancestors` built from the same map. Unknown terms yield the
    empty set, which the caller reads as "nothing to test against".
    """

    def parents(term: str) -> Set[str]:
        return set(child_to_parents.get(term, ()))

    return parents


def propagate_annotation_map(
    protein_terms: Dict[str, Set[str]],
    ancestors_fn: Callable[[str], Iterable[str]],
) -> Dict[str, Set[str]]:
    """Apply the True Path Rule to a ``{protein: {terms}}`` annotation map.

    Each protein's term set is closed over ``ancestors_fn``: an annotation to a
    child term implies its ancestors by definition. Ancestors are looked up
    once per distinct term, not once per (protein, term) — the annotation map
    has ~10^6 pairs over ~10^4 terms. The input map is not mutated.
    """
    cache: Dict[str, frozenset] = {}
    propagated: Dict[str, Set[str]] = {}
    for protein, terms in protein_terms.items():
        closure = set(terms)
        for term in terms:
            cached = cache.get(term)
            if cached is None:
                cached = cache[term] = frozenset(ancestors_fn(term))
            closure |= cached
        propagated[protein] = closure
    return propagated


def nearest_parents(
    ordered_ancestors_fn: Callable[[str], Sequence[str]],
) -> Callable[[str], Set[str]]:
    """Direct-parents function for a hierarchy encoded in the term id.

    :func:`dotted_ancestors`, :func:`alpha_prefix_ancestors` and
    :func:`src.ec_annotation_source.ec_ancestors` all return ancestors **most
    specific first**, and these classifications are trees — every id has exactly
    one parent — so the direct parent is simply the first entry.
    """

    def parents(term: str) -> Set[str]:
        ancestors = ordered_ancestors_fn(term)
        return {ancestors[0]} if ancestors else set()

    return parents


def propagate_via_ancestors(
    direct_associations: Iterable[Any],
    ancestors_fn: Callable[[str], Iterable[str]],
) -> "List[Annotation]":
    """Propagate associations with independent provenance and evidence merging.

    Every direct association supports its own ``(domain, term)`` pair and every
    ancestor returned by ``ancestors_fn``. Supporting evidence is aggregated
    independently of provenance:

    * a pair is labelled ``direct`` whenever direct evidence exists, and its
      ``direct_source_term`` is the pair's own term;
    * otherwise the propagated source with the lowest q-value supplies
      ``direct_source_term`` (term id breaks equal-q ties);
    * ``q_value`` is the minimum and ``association_score`` the maximum across
      all direct and propagated support for the pair.

    Thus a strong child may strengthen a directly observed parent without
    relabelling that parent as propagated. Results do not depend on input order.

    Args:
        direct_associations: iterable of objects exposing ``domain``,
            ``go_term``, ``q_value`` and ``hyper_score``.
        ancestors_fn: maps a term id to its ancestor term ids.

    Returns:
        De-duplicated ``ontology_processor.Annotation`` objects.
    """
    from src.ontology_processor import Annotation

    ordered = sorted(
        direct_associations,
        key=lambda assoc: (
            assoc.q_value,
            assoc.domain,
            assoc.go_term,
            -assoc.hyper_score,
        ),
    )

    # Each support record is (association, is_direct). Preserve a stable pair
    # order: direct pairs first, followed by newly encountered ancestors.
    support: Dict[tuple[str, str], List[tuple[Any, bool]]] = defaultdict(list)
    pair_order: List[tuple[str, str]] = []

    def add_support(key: tuple[str, str], assoc: Any, is_direct: bool) -> None:
        if key not in support:
            pair_order.append(key)
        support[key].append((assoc, is_direct))

    for assoc in ordered:
        add_support((assoc.domain, assoc.go_term), assoc, True)

    for assoc in ordered:
        for ancestor in sorted(ancestors_fn(assoc.go_term)):
            add_support((assoc.domain, ancestor), assoc, False)

    annotations: List[Annotation] = []
    for domain, term in pair_order:
        evidence = support[(domain, term)]
        direct = [assoc for assoc, is_direct in evidence if is_direct]
        if direct:
            annotation_type = "direct"
            direct_source_term = term
        else:
            source = min(
                (assoc for assoc, _ in evidence),
                key=lambda assoc: (assoc.q_value, assoc.go_term),
            )
            annotation_type = "propagated"
            direct_source_term = source.go_term

        annotations.append(
            Annotation(
                domain=domain,
                go_term=term,
                q_value=min(assoc.q_value for assoc, _ in evidence),
                association_score=max(assoc.hyper_score for assoc, _ in evidence),
                annotation_type=annotation_type,
                direct_source_term=direct_source_term,
            )
        )

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
