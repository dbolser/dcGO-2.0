"""Relative inference: the parental-background enrichment test, ontology-agnostic.

This is Step 2's second statistical inference in the dcGO paper (Fang & Gough
2013) — *not* the True Path Rule, which is Step 3 and lives in
:mod:`src.hierarchy`. The overall inference asks whether a domain is enriched
for a term against the whole analysable background; this one asks whether it is
*still* enriched against only those proteins annotated to the term's **direct
parents**.

Both questions have to be answered, because each alone fails in a
characteristic way:

* Overall alone cannot locate the level of the hierarchy the signal lives at. A
  domain genuinely associated with a broad term makes every descendant look
  enriched for free, since the descendant's proteins are a subset of the
  parent's.
* Relative alone rests on a background that can be small and idiosyncratic, and
  root terms have no parents to test against at all.

The module implements the *test*. How its verdict is combined with the overall
p-value is the caller's business, and the two current callers differ:

* ``run_dcgo_human.py --enable-relative-inference`` applies it as a post-hoc
  ``alpha`` filter after BH. **This is not the paper's method** — the paper
  takes ``max(overall_p, relative_p)`` (an intersection-union statistic, valid
  as a p-value without correction) and applies BH to that. Tracked as
  ``VALIDATION_PLAN.md`` next-steps item 2.
* Callers wanting paper parity should use :func:`relative_p_values` and combine
  before correcting.

Hierarchy access is injected as two functions — ``parents_fn`` (direct parents,
what the test ranges over) and ``ancestors_fn`` (transitive, what the background
index propagates over) — so every ontology in the registry with a hierarchy can
use this, not only GO.
"""

from __future__ import annotations

from typing import Callable, Dict, Iterable, List, Sequence, Set, Tuple

from loguru import logger
from scipy.stats import fisher_exact

ParentsFn = Callable[[str], Iterable[str]]
AncestorsFn = Callable[[str], Iterable[str]]


class InsufficientBackgroundError(ValueError):
    """The parental statistical universe is smaller than the configured minimum."""


class InvalidContingencyTableError(ValueError):
    """The parental-background counts cannot form a valid contingency table."""


class BackgroundIndex:
    """Inverted term→proteins and domain→proteins indexes for the parental test.

    Term membership is **propagated**: a protein annotated to a term is counted
    under every ancestor of that term. That is the True Path Rule, and it is
    what "the proteins annotated to the parent" has to mean for the test to be
    answerable — indexing only direct annotations gives any parent nobody is
    directly annotated to an empty background, so the test raises and every
    child of it is discarded untested. Domain membership is direct, because
    domains have no hierarchy here.

    Inverting the maps is also what makes the stage finish at all: computed
    straight from the maps, each cell is a pass over the whole proteome *per
    parent test*.

    See ``docs/design/ontology-processing.md`` for the regression history.
    """

    __slots__ = ("term_proteins", "domain_proteins")

    def __init__(
        self,
        protein_domain_map: Dict[str, List[str]],
        protein_terms: Dict[str, Set[str]],
        ancestors_fn: AncestorsFn | None = None,
    ) -> None:
        term_proteins: Dict[str, set] = {}
        # Ancestors are looked up once per distinct term, not once per
        # (protein, term): the annotation map has ~10^6 pairs over ~10^4 terms.
        ancestor_cache: Dict[str, Set[str]] = {}
        for protein, terms in protein_terms.items():
            for term in terms:
                term_proteins.setdefault(term, set()).add(protein)
                if ancestors_fn is None:
                    continue
                if term not in ancestor_cache:
                    ancestor_cache[term] = set(ancestors_fn(term))
                for ancestor in ancestor_cache[term]:
                    term_proteins.setdefault(ancestor, set()).add(protein)

        domain_proteins: Dict[str, set] = {}
        for protein, domains in protein_domain_map.items():
            for domain in domains:
                domain_proteins.setdefault(domain, set()).add(protein)

        self.term_proteins = term_proteins
        self.domain_proteins = domain_proteins


def parental_p_value(
    domain: str,
    child_term: str,
    parent_term: str,
    index: BackgroundIndex,
    min_background_size: int,
) -> float:
    """One-tailed Fisher p for *domain* × *child_term* within *parent_term*.

    Raises:
        InsufficientBackgroundError: the parent background is below the minimum.
        InvalidContingencyTableError: the counts cannot form a 2x2 table.
    """
    parent_proteins = index.term_proteins.get(parent_term, frozenset())
    if len(parent_proteins) < min_background_size:
        raise InsufficientBackgroundError(
            f"Insufficient background size: {len(parent_proteins)} "
            f"< {min_background_size}"
        )

    # Set algebra on the index. Identical counts to scanning every protein.
    domain_proteins = index.domain_proteins.get(domain, frozenset())
    child_in_background = index.term_proteins.get(child_term, frozenset()) & (
        parent_proteins
    )
    domain_in_background = domain_proteins & parent_proteins

    a = len(child_in_background & domain_proteins)
    b = len(child_in_background) - a
    c = len(domain_in_background) - a
    d = len(parent_proteins) - (a + b + c)

    if a == 0 or (a + b) == 0 or (a + c) == 0:
        return 1.0  # no association possible

    if d < 0:
        raise InvalidContingencyTableError(
            f"Invalid contingency table: a={a}, b={b}, c={c}, d={d}"
        )

    try:
        _, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        return float(p_value)
    except (ValueError, ZeroDivisionError) as exc:
        logger.warning(f"Fisher's exact test failed: {exc}")
        return 1.0


def relative_p_value(
    domain: str,
    term: str,
    index: BackgroundIndex,
    parents_fn: ParentsFn,
    min_background_size: int,
    reject_at: float | None = None,
) -> float:
    """The association's relative p-value: its *weakest* parental result.

    The association must survive *every* direct parent, so the governing
    p-value is the largest — the same intersection-union logic the paper applies
    when combining the overall and relative inferences.

    A term with no parents (a root, or one absent from the hierarchy) has no
    parental background to test against and returns 0.0, which is the identity
    for the caller's ``max(overall_p, relative_p)`` combination — an untestable
    term is decided by the overall inference alone rather than being discarded.

    Args:
        reject_at: if given, stop at the first parent scoring at or above this
            and return that p-value. The filter decision is unchanged (any
            parent at or above the threshold puts the maximum there too), but it
            skips the remaining parents — which preserves both the cost and the
            rejection tally of the original short-circuiting implementation.
            Leave ``None`` to always compute the true maximum, which is what a
            caller combining before BH needs.

    Raises:
        InsufficientBackgroundError, InvalidContingencyTableError: propagated
            from :func:`parental_p_value` for the caller's rejection policy.
    """
    worst = 0.0
    for parent in parents_fn(term):
        p_value = parental_p_value(domain, term, parent, index, min_background_size)
        worst = max(worst, p_value)
        if reject_at is not None and p_value >= reject_at:
            return p_value
    return worst


def relative_p_values(
    associations: Sequence,
    protein_domain_map: Dict[str, List[str]],
    protein_terms: Dict[str, Set[str]],
    parents_fn: ParentsFn,
    ancestors_fn: AncestorsFn,
    min_background_size: int = 3,
    reject_at: float | None = None,
) -> Tuple[List[float | None], Dict[str, int]]:
    """Relative p-value per association, ``None`` where the test could not run.

    Returns ``(p_values, rejection_counts)``. ``rejection_counts`` is keyed by
    exception type name, for the caller to report in aggregate — at ~10^5
    occurrences, logging each one produced a 19 MB log file.

    ``reject_at`` is forwarded to :func:`relative_p_value`; pass it only when
    the result feeds a threshold decision, never when it feeds a combination.
    """
    index = BackgroundIndex(protein_domain_map, protein_terms, ancestors_fn)
    p_values: List[float | None] = []
    rejections: Dict[str, int] = {}

    for i, assoc in enumerate(associations):
        if i and i % 10000 == 0:
            logger.info(f"  relative inference {i:,}/{len(associations):,}")
        try:
            p_values.append(
                relative_p_value(
                    assoc.domain,
                    assoc.go_term,
                    index,
                    parents_fn,
                    min_background_size,
                    reject_at=reject_at,
                )
            )
        except (InsufficientBackgroundError, InvalidContingencyTableError) as exc:
            # Data-dependent, expected, and counted. Programming errors are not
            # caught here and must remain visible.
            name = type(exc).__name__
            rejections[name] = rejections.get(name, 0) + 1
            p_values.append(None)

    return p_values, rejections


def filter_by_parental_background(
    associations: Sequence,
    protein_domain_map: Dict[str, List[str]],
    protein_terms: Dict[str, Set[str]],
    parents_fn: ParentsFn,
    ancestors_fn: AncestorsFn,
    min_background_size: int = 3,
    alpha_threshold: float = 0.05,
) -> List:
    """Keep only associations still enriched within their parents' background.

    The post-hoc form of the test: associations whose relative p-value exceeds
    ``alpha_threshold`` are dropped, as are those the test could not evaluate
    (conservative — an association that cannot be shown specific is not kept).

    Prefer :func:`relative_p_values` plus ``max(overall_p, relative_p)`` before
    BH for paper parity; see the module docstring.
    """
    if not associations:
        logger.warning("No associations provided for relative inference")
        return []
    if not protein_domain_map or not protein_terms:
        raise ValueError("Protein mapping data cannot be empty")

    p_values, rejections = relative_p_values(
        associations,
        protein_domain_map,
        protein_terms,
        parents_fn,
        ancestors_fn,
        min_background_size,
        reject_at=alpha_threshold,
    )

    # Strictly below the threshold, matching the original ``p >= alpha`` reject.
    kept = [
        assoc
        for assoc, p in zip(associations, p_values)
        if p is not None and p < alpha_threshold
    ]

    logger.info(
        f"Relative inference: {len(kept)}/{len(associations)} associations retained"
    )
    if rejections:
        detail = ", ".join(f"{n:,} x {kind}" for kind, n in sorted(rejections.items()))
        logger.info(
            f"  {sum(rejections.values()):,} parent tests could not be evaluated "
            f"and their associations were rejected untested ({detail}). The "
            "background is propagated, so this should be rare; a large count "
            "means the hierarchy or the annotation map is not what the test "
            "expects."
        )
    return kept
