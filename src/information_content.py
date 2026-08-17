"""Annotation-frequency information content of ontology terms.

One IC convention for the whole project: ``IC(t) = -log2(P(t))`` where ``P(t)``
is the fraction of proteins in the analysable universe annotated to ``t``. It
is the marginal (CAFA-style) IC the temporal benchmark has always used
(``validation/temporal_benchmark.py``), now shared so the pipeline's reported
IC and the benchmark's IC cells cannot drift apart.

Two properties matter to callers:

* **The frequencies must come from a True-Path-propagated map** whenever the
  ontology has a hierarchy. An annotation to a child term implies its
  ancestors, so an unpropagated ``P(t)`` understates the frequency of every
  non-leaf term and inflates its IC. For a flat vocabulary the direct map *is*
  the propagated map. The caller owns propagation (see
  :func:`src.hierarchy.propagate_annotation_map`); this module just counts.
* **Roots have IC 0 by construction**: a term carried by every protein has
  ``P(t) = 1``. That is what makes IC a usable floor against vacuous
  associations to the top of a DAG — any positive floor removes them.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, Iterable, Mapping


def information_content(
    annotation_map: Mapping[str, Iterable[str]],
) -> Dict[str, float]:
    """Marginal IC per term from a ``{protein: {terms}}`` annotation map.

    ``P(t)`` is estimated over the proteins in *annotation_map* — the run's own
    analysable universe — and ``IC(t) = -log2(P(t))``. Terms with ``P(t) = 1``
    (roots, universal terms) get IC 0.0 exactly. Terms absent from the map are
    absent from the result: they carry no frequency estimate, and callers
    should treat them as IC 0 ("no information"), matching the temporal
    benchmark's convention.

    The map must already be True-Path propagated when the ontology has a
    hierarchy — see the module docstring.
    """
    n_proteins = len(annotation_map)
    if n_proteins == 0:
        return {}
    counts: Dict[str, int] = defaultdict(int)
    for terms in annotation_map.values():
        for term in terms:
            counts[term] += 1
    ic: Dict[str, float] = {}
    for term, count in counts.items():
        p = count / n_proteins
        ic[term] = -math.log2(p) if 0.0 < p < 1.0 else 0.0
    return ic
