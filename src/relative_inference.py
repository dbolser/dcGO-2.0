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

The paper combines the two inferences by taking ``max(overall_p, relative_p)``
and applying BH to *that*. The maximum is an intersection-union statistic: the
null is "fails at least one inference", rejecting requires rejecting both, and
the maximum of the individual p-values is itself a valid p-value for that null
with no multiplicity correction (Berger 1982). So it is exactly the quantity BH
is entitled to correct — which is why the combination has to happen *before* the
correction, not as a filter afterwards.

Two implementations of the same test live here:

* :func:`compute_relative_p_values` — vectorised over every co-occurring pair,
  via sparse matmuls. This is what ``run_dcgo_human.py`` uses, because paper
  parity needs a relative p-value for all ~10^6-10^7 candidate pairs.
* :func:`relative_p_value` and friends — the scalar, set-algebra reference
  implementation. Kept because it is obviously correct by inspection and the
  vectorised path is validated against it, pair for pair, in the tests.

:func:`filter_by_parental_background` applies the test as a post-hoc ``alpha``
filter instead. That is **not** the paper's method and is no longer reachable
from the CLI; it survives only as the scalar path's driver.

Hierarchy access is injected as two functions — ``parents_fn`` (direct parents,
what the test ranges over) and ``ancestors_fn`` (transitive, what the background
index propagates over) — so every ontology in the registry with a hierarchy can
use this, not only GO.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable, Dict, Iterable, List, Sequence, Set, Tuple

import numpy as np

if TYPE_CHECKING:
    from scipy import sparse

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


# ---------------------------------------------------------------------------- #
# Vectorised path: relative p-values for every enumerated pair, before BH       #
# ---------------------------------------------------------------------------- #
#
# The loop above tests a handful of already-significant associations. Paper
# parity needs the relative p-value for *every* co-occurring pair, because the
# quantity BH corrects is ``max(overall_p, relative_p)`` — so the loop's ~10^4
# scipy calls become ~10^6-10^7 and have to be vectorised.
#
# The algebra collapses nicely. Under the True Path Rule a child's proteins are
# a subset of its parent's, so within the parent background:
#
#     a = |proteins(d) & proteins(c)|            (propagated co-occurrence)
#     b = |proteins(c)| - a
#     c = |proteins(d) & proteins(p)| - a        (propagated co-occurrence)
#     d = |proteins(p)| - |proteins(c)| - |proteins(d) & proteins(p)| + a
#
# Every term is either a propagated column sum or an entry of the propagated
# domain x term co-occurrence product — the same sparse matmul the overall test
# already uses, just against a propagated term matrix.


def build_propagated_term_matrix(
    protein_term_matrix: "sparse.csr_matrix",
    term_ids: Sequence[str],
    ancestors_fn: AncestorsFn,
) -> Tuple["sparse.csr_matrix", List[str], Dict[str, int]]:
    """Propagate a protein x term matrix up the hierarchy.

    The term axis is *extended* with every ancestor of an annotated term, even
    ancestors nothing is annotated to directly. That is the whole point: a
    parent with no direct annotation still needs a background, and building it
    from the unpropagated map is the defect that made the filter reject 54,951
    associations untested (#46).

    The extension is internal to this test. The hypothesis space BH corrects
    over stays exactly ``term_ids``.

    Returns:
        ``(propagated, extended_ids, index_of)`` — a binary
        (n_proteins, n_extended) matrix, the extended term list, and its
        term → column lookup.
    """
    from scipy import sparse

    extended: List[str] = list(term_ids)
    index_of: Dict[str, int] = {term: i for i, term in enumerate(extended)}
    for term in term_ids:
        for ancestor in ancestors_fn(term):
            if ancestor not in index_of:
                index_of[ancestor] = len(extended)
                extended.append(ancestor)

    # Incidence matrix: row t marks t itself and all of its ancestors.
    rows: List[int] = []
    cols: List[int] = []
    for i, term in enumerate(term_ids):
        rows.append(i)
        cols.append(i)
        for ancestor in ancestors_fn(term):
            rows.append(i)
            cols.append(index_of[ancestor])
    incidence = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.int32), (rows, cols)),
        shape=(len(term_ids), len(extended)),
    )

    propagated = (protein_term_matrix.astype(np.int32) @ incidence).tocsr()
    # A protein reaching one ancestor by several routes must still count once.
    propagated.data[:] = 1
    return propagated, extended, index_of


def compute_relative_p_values(
    protein_domain_matrix: "sparse.csr_matrix",
    protein_term_matrix: "sparse.csr_matrix",
    pair_index: "np.ndarray",
    term_ids: Sequence[str],
    parents_fn: ParentsFn,
    ancestors_fn: AncestorsFn,
    min_background_size: int = 3,
    domain_block: int = 4096,
    fisher_batch_size: int = 50000,
) -> Tuple["np.ndarray", "np.ndarray", Dict[str, int]]:
    """Relative p-value and governing 2x2 table for every enumerated pair.

    Args:
        pair_index: dense domain-major indices as returned by
            ``compute_cooccurring_contingency_tables``.

    Returns:
        ``(relative_p, tables, rejections)``. ``relative_p[i]`` is the largest
        p-value over the pair's direct parents — the association must survive
        every one of them, the same intersection-union logic the paper uses to
        combine the overall and relative inferences. A pair whose term has no
        parents scores 0.0, so ``max(overall_p, relative_p)`` leaves it to the
        overall inference. A pair with any untestable parent scores 1.0, which
        is the vectorised spelling of the loop's "conservatively reject".
        ``tables[i]`` is the governing parent's table (zeros where untested).
    """
    from scipy import sparse

    from src.vectorized_fisher import fisher_exact_parallel

    n_terms = len(term_ids)
    n_pairs = len(pair_index)

    propagated, extended_ids, index_of = build_propagated_term_matrix(
        protein_term_matrix, term_ids, ancestors_fn
    )
    logger.info(
        f"  propagated term axis: {n_terms:,} annotated -> "
        f"{len(extended_ids):,} with ancestors"
    )
    term_counts = np.asarray(propagated.sum(axis=0)).ravel().astype(np.int64)

    # Propagated domain x term co-occurrence, blocked over domains for the same
    # allocation reason compute_cooccurring_contingency_tables blocks.
    domain_csc = protein_domain_matrix.astype(np.int32).tocsc()
    n_domains = protein_domain_matrix.shape[1]
    blocks = []
    for start in range(0, n_domains, domain_block):
        stop = min(start + domain_block, n_domains)
        blocks.append((domain_csc[:, start:stop].T @ propagated).tocsr())
    cooccurrence = sparse.vstack(blocks, format="csr") if blocks else None
    del blocks, domain_csc

    domain_idx = pair_index // n_terms
    term_idx = pair_index % n_terms

    # Flatten (pair, parent) into one axis, remembering each parent's pair.
    parent_cols: List[int] = []
    owner: List[int] = []
    for i in range(n_pairs):
        for parent in parents_fn(term_ids[term_idx[i]]):
            column = index_of.get(parent)
            if column is not None:
                parent_cols.append(column)
                owner.append(i)

    relative_p = np.zeros(n_pairs, dtype=np.float64)  # no parents -> 0.0
    tables = np.zeros((n_pairs, 2, 2), dtype=np.int32)
    rejections: Dict[str, int] = {}
    if not owner:
        return relative_p, tables, rejections

    owner_arr = np.asarray(owner, dtype=np.int64)
    parent_arr = np.asarray(parent_cols, dtype=np.int64)
    del owner, parent_cols

    child_arr = term_idx[owner_arr]
    dom_arr = domain_idx[owner_arr]

    a = np.asarray(cooccurrence[dom_arr, child_arr]).ravel().astype(np.int64)
    dp = np.asarray(cooccurrence[dom_arr, parent_arr]).ravel().astype(np.int64)
    n_child = term_counts[child_arr]
    n_parent = term_counts[parent_arr]

    b = n_child - a
    c = dp - a
    d = n_parent - n_child - dp + a

    # Guards, in the same order and with the same verdicts as the loop.
    too_small = n_parent < min_background_size
    invalid = d < 0
    untestable = too_small | invalid
    if too_small.any():
        rejections["InsufficientBackgroundError"] = int(too_small.sum())
    if invalid.any():
        rejections["InvalidContingencyTableError"] = int(invalid.sum())

    # a == 0 or an empty margin means no association is possible: p = 1 exactly.
    trivial = (a == 0) | ((a + b) == 0) | ((a + c) == 0)

    p_per_parent = np.ones(len(owner_arr), dtype=np.float64)
    testable = ~(untestable | trivial)
    if testable.any():
        sub = np.empty((int(testable.sum()), 2, 2), dtype=np.int32)
        sub[:, 0, 0] = a[testable]
        sub[:, 0, 1] = b[testable]
        sub[:, 1, 0] = c[testable]
        sub[:, 1, 1] = d[testable]
        _, p_sub = fisher_exact_parallel(
            sub, alternative="greater", batch_size=fisher_batch_size
        )
        p_per_parent[testable] = p_sub

    # An untestable parent is conservatively fatal for its association.
    p_per_parent[untestable] = 1.0

    # Group maximum, and the governing parent's table alongside it.
    order = np.lexsort((p_per_parent, owner_arr))
    sorted_owner = owner_arr[order]
    last = np.ones(len(order), dtype=bool)
    last[:-1] = sorted_owner[:-1] != sorted_owner[1:]
    winners = order[last]
    winning_pairs = sorted_owner[last]

    relative_p[winning_pairs] = p_per_parent[winners]
    tables[winning_pairs, 0, 0] = a[winners]
    tables[winning_pairs, 0, 1] = b[winners]
    tables[winning_pairs, 1, 0] = c[winners]
    tables[winning_pairs, 1, 1] = np.maximum(d[winners], 0)

    return relative_p, tables, rejections
