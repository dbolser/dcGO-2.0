"""Calibrated knowledge projection through shared units.

dcGO's core idea, unchanged: function is carried by reusable units — domains,
and the combinations they occur in — so what is known about the proteins that
carry a unit can be projected onto proteins nobody has studied. This module
changes *what is projected*.

dcGO projects a significance: a Fisher p-value says an association is
unlikely to be chance, not how often it holds. A domain seen in 1,000
annotated proteins, 100 of which have a term, is overwhelmingly significant
and right one time in ten; a domain seen in 4 proteins, all with the term, is
barely significant and nearly always right. The p-score then min-max scales
each protein's scores to [0, 1], so a "1.0" on one protein and a "1.0" on
another carry no shared meaning and the predictions cannot be pooled into one
list a curator could work down.

Here every (unit, term) pair carries an estimate of the quantity projection
needs, ``P(term | protein carries unit)``, estimated on the proteins
annotated **in that aspect** (a protein with no biological-process
annotation is not evidence against any biological-process term) and shrunk
towards the term's base rate with a Beta prior of strength ``s``::

    θ(u, t) = (k(u, t) + s · π(t)) / (n(u) + s)

``n(u)`` carriers of ``u`` are annotated in the aspect and ``k(u, t)`` of them
have ``t`` (True-Path propagated). A protein's score for ``t`` is the best
``θ`` over its units: the most specific unit it shares with annotated
proteins speaks for it, and overlapping InterPro entries for one region
cannot add up. With propagated labels and ``π(parent) ≥ π(child)`` the scores
are hierarchy-consistent by construction (a parent never scores below its
child), and they are on one probability scale across proteins.

The same estimator takes any unit family: InterPro entries, the bag-of-domain
pairs and triples dcGO calls supra-domains, a protein's whole distinct-domain
architecture, member-database signatures (PANTHER subfamilies included).
:func:`knn_project` is the literature counterpart: proteins are near each
other when the papers about them are, and a protein inherits its literature
neighbours' terms.
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import scipy.sparse as sp

#: GO aspect roots; never predicted or scored (CAFA convention).
ASPECT_ROOTS = frozenset({"GO:0008150", "GO:0003674", "GO:0005575"})


# --------------------------------------------------------------------------- #
# Labels                                                                       #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class AspectTraining:
    """Propagated training labels for one GO aspect.

    ``proteins`` are exactly the proteins with at least one non-root term in
    the aspect; only they say anything about the aspect's terms.
    """

    aspect: str
    proteins: tuple[str, ...]
    terms: tuple[str, ...]
    labels: sp.csr_matrix  # proteins x terms, 1.0 where annotated
    prior: np.ndarray  # fraction of ``proteins`` carrying each term

    @property
    def term_index(self) -> dict[str, int]:
        return {t: i for i, t in enumerate(self.terms)}


def build_aspect_training(
    annotation_map: Mapping[str, Iterable[str]],
    get_ancestors,
    term_aspect: Mapping[str, str],
    aspect: str,
    exclude: Iterable[str] = (),
) -> AspectTraining:
    """Labels for ``aspect`` from an unpropagated ``{protein: terms}`` map."""
    excluded = set(exclude)
    closures: dict[str, set[str]] = {}
    for protein, terms in annotation_map.items():
        if protein in excluded:
            continue
        closed: set[str] = set()
        for t in terms:
            closed.add(t)
            closed.update(get_ancestors(t))
        kept = {t for t in closed if term_aspect.get(t) == aspect} - ASPECT_ROOTS
        if kept:
            closures[protein] = kept
    proteins = tuple(sorted(closures))
    terms = tuple(sorted({t for ts in closures.values() for t in ts}))
    t_index = {t: i for i, t in enumerate(terms)}
    rows, cols = [], []
    for i, p in enumerate(proteins):
        for t in closures[p]:
            rows.append(i)
            cols.append(t_index[t])
    labels = sp.csr_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(len(proteins), len(terms))
    )
    prior = np.asarray(labels.sum(axis=0)).ravel() / max(len(proteins), 1)
    return AspectTraining(aspect, proteins, terms, labels, prior)


def build_aspect_trainings(
    annotation_map: Mapping[str, Iterable[str]],
    get_ancestors,
    term_aspect: Mapping[str, str],
    aspects: Sequence[str],
    exclude: Mapping[str, Iterable[str]] | None = None,
) -> dict[str, AspectTraining]:
    """:func:`build_aspect_training` for several aspects in one sparse product.

    The True-Path closure is ``direct @ ancestor-or-self`` over every term
    used, which stays fast for an all-species map of millions of proteins.
    ``exclude[aspect]`` lists proteins kept out of that aspect's training.
    Same result as calling :func:`build_aspect_training` per aspect.
    """
    proteins = sorted(annotation_map)
    used = sorted({t for ts in annotation_map.values() for t in ts})
    closure_terms = sorted(set(used).union(*(get_ancestors(t) for t in used)))
    t_index = {t: i for i, t in enumerate(closure_terms)}
    rows, cols = [], []
    for i, p in enumerate(proteins):
        for t in annotation_map[p]:
            rows.append(i)
            cols.append(t_index[t])
    direct = sp.csr_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(len(proteins), len(closure_terms))
    )
    closed = (direct @ ancestor_matrix(closure_terms, get_ancestors)).tocsc()
    closed.data[:] = 1.0
    out: dict[str, AspectTraining] = {}
    for aspect in aspects:
        cols_a = [
            j
            for j, t in enumerate(closure_terms)
            if term_aspect.get(t) == aspect and t not in ASPECT_ROOTS
        ]
        sub = closed[:, cols_a].tocsr()
        keep = np.diff(sub.indptr) > 0
        if exclude and aspect in exclude:
            excluded = set(exclude[aspect])
            keep &= np.array([p not in excluded for p in proteins])
        sub = sub[keep]
        present = np.asarray(sub.sum(axis=0)).ravel() > 0
        sub = sub[:, present].tocsr()
        terms = tuple(np.asarray([closure_terms[j] for j in cols_a])[present])
        kept = tuple(np.asarray(proteins, dtype=object)[keep])
        prior = np.asarray(sub.sum(axis=0)).ravel() / max(len(kept), 1)
        out[aspect] = AspectTraining(aspect, kept, terms, sub, prior)
    return out


# --------------------------------------------------------------------------- #
# Units                                                                        #
# --------------------------------------------------------------------------- #
class UnitVocabulary:
    """Index of unit ids, and protein x unit incidence matrices over it."""

    def __init__(self, protein_units: Mapping[str, Iterable[str]]):
        self.protein_units = {
            p: tuple(dict.fromkeys(u)) for p, u in protein_units.items()
        }
        units = sorted({u for us in self.protein_units.values() for u in us})
        self.units: tuple[str, ...] = tuple(units)
        self.index = {u: i for i, u in enumerate(units)}

    def incidence(self, proteins: Sequence[str]) -> sp.csr_matrix:
        rows, cols = [], []
        for i, p in enumerate(proteins):
            for u in self.protein_units.get(p, ()):
                rows.append(i)
                cols.append(self.index[u])
        return sp.csr_matrix(
            (np.ones(len(rows)), (rows, cols)), shape=(len(proteins), len(self.units))
        )


def read_protein_signatures(
    protein2ipr: Path, proteins: set[str] | None = None, prefix: str = "SIG:"
) -> dict[str, set[str]]:
    """Member-database signatures per protein (column 4 of ``protein2ipr``).

    PANTHER subfamilies (``PTHR…:SF…``) are kept as their own units, and
    their family is added too, so both levels are available.
    """
    out: dict[str, set[str]] = defaultdict(set)
    with gzip.open(protein2ipr, "rt") as handle:
        for line in handle:
            fields = line.split("\t", 4)
            if len(fields) < 4:
                continue
            acc, sig = fields[0], fields[3]
            if proteins is not None and acc not in proteins:
                continue
            out[acc].add(prefix + sig)
            if ":SF" in sig:
                out[acc].add(prefix + sig.split(":", 1)[0])
    return dict(out)


def architecture_units(
    protein_domains: Mapping[str, Iterable[str]], min_size: int = 4
) -> dict[str, set[str]]:
    """One unit per protein: its whole distinct-domain set, when it has ``≥ min_size``.

    Pairs and triples are already supra-domain units; this adds the exact
    architecture of larger proteins, the most specific context dcGO's
    combinations stop short of. Only single-domain ids (no ``,``) are used.
    """
    out: dict[str, set[str]] = {}
    for p, domains in protein_domains.items():
        singles = sorted({d for d in domains if "," not in d})
        if len(singles) >= min_size:
            out[p] = {"ARCH:" + "+".join(singles)}
    return out


def bag_supra_units(
    protein_domains: Mapping[str, Iterable[str]], max_size: int = 3
) -> dict[str, set[str]]:
    """All unordered pairs/triples of a protein's distinct domains (``bag`` mode)."""
    out: dict[str, set[str]] = {}
    for p, domains in protein_domains.items():
        singles = sorted({d for d in domains if "," not in d})
        units: set[str] = set()
        for k in range(2, max_size + 1):
            units.update(",".join(c) for c in combinations(singles, k))
        if units:
            out[p] = units
    return out


def unit_family(unit: str) -> str:
    """Which family a unit id belongs to, from its spelling."""
    if unit.startswith("ARCH:"):
        return "arch"
    if unit.startswith("SIG:"):
        return "sig"
    return "supra" if "," in unit else "ipr"


def split_families(
    units: Mapping[str, Iterable[str]],
) -> dict[str, dict[str, set[str]]]:
    """``{family: {protein: units}}`` from a mixed ``{protein: units}`` map."""
    out: dict[str, dict[str, set[str]]] = {
        f: {} for f in ("ipr", "supra", "arch", "sig")
    }
    for p, us in units.items():
        for u in us:
            out[unit_family(u)].setdefault(p, set()).add(u)
    return out


def read_restricted_units(
    protein2ipr: Path,
    proteins: set[str],
    vocabulary: set[str],
    min_domain_length: int = 10,
) -> dict[str, set[str]]:
    """Units of ``proteins`` that are in ``vocabulary``, straight from ``protein2ipr``.

    For training on a large background (all species): only units some query
    protein carries can ever be projected, so only those are kept, and
    supra-domains are generated only from a protein's in-vocabulary entries.
    Entries follow :class:`DomainAnnotationParser` (rows shorter than
    ``min_domain_length`` dropped; supra-domains are sorted bag combinations
    of 2–3 distinct entries); signatures follow
    :func:`read_protein_signatures` (every row).
    """
    singles: dict[str, set[str]] = defaultdict(set)
    sigs: dict[str, set[str]] = defaultdict(set)
    constituents = {d for u in vocabulary if "," in u for d in u.split(",")}
    with gzip.open(protein2ipr, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[0] not in proteins:
                continue
            acc, ipr, sig = fields[0], fields[1], fields[3]
            for unit in ("SIG:" + sig, "SIG:" + sig.split(":", 1)[0]):
                if unit in vocabulary:
                    sigs[acc].add(unit)
            if int(fields[5]) - int(fields[4]) + 1 >= min_domain_length:
                singles[acc].add(ipr)
    out: dict[str, set[str]] = {}
    for acc in set(singles) | set(sigs):
        doms = singles.get(acc, set())
        units = {d for d in doms if d in vocabulary} | sigs.get(acc, set())
        relevant = sorted(d for d in doms if d in constituents)
        for k in (2, 3):
            units.update(
                c
                for c in (",".join(x) for x in combinations(relevant, k))
                if c in vocabulary
            )
        if len(doms) >= 4:
            arch = "ARCH:" + "+".join(sorted(doms))
            if arch in vocabulary:
                units.add(arch)
        if units:
            out[acc] = units
    return out


def merge_units(*families: Mapping[str, Iterable[str]]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for family in families:
        for p, units in family.items():
            out[p].update(units)
    return dict(out)


# --------------------------------------------------------------------------- #
# Unit -> term evidence                                                        #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class UnitTermEvidence:
    """Counts behind ``θ(u, t)`` for one aspect."""

    n_carriers: np.ndarray  # per unit: annotated-in-aspect carriers
    counts: sp.csr_matrix  # unit x term: carriers with the term
    prior: np.ndarray  # per term

    def posterior(
        self, strength: float, quantile: float | None = None
    ) -> sp.csr_matrix:
        """``θ`` on the non-zero counts (pairs never observed carry no evidence).

        By default the posterior mean. With ``quantile`` (e.g. 0.1) the lower
        credible bound of ``Beta(k + s·π, n − k + s·(1 − π))`` instead: a unit
        seen in one or two annotated proteins then cannot outbid a unit with
        the same rate over hundreds — the robust choice when the best unit
        speaks for the protein.
        """
        c = self.counts.tocoo()
        a = c.data + strength * self.prior[c.col]
        b = self.n_carriers[c.row] - c.data + strength * (1.0 - self.prior[c.col])
        if quantile is None:
            theta = a / (a + b)
        else:
            from scipy.special import betaincinv

            theta = betaincinv(a, np.maximum(b, 1e-9), quantile)
        return sp.csr_matrix((theta, (c.row, c.col)), shape=c.shape)


def unit_term_evidence(
    train_incidence: sp.csr_matrix, training: AspectTraining
) -> UnitTermEvidence:
    """``n(u)`` and ``k(u, t)`` from the training proteins' units and labels."""
    if train_incidence.shape[0] != training.labels.shape[0]:
        raise ValueError("incidence rows must be the training proteins, in order")
    n = np.asarray(train_incidence.sum(axis=0)).ravel()
    counts = (train_incidence.T @ training.labels).tocsr()
    counts.eliminate_zeros()
    return UnitTermEvidence(n_carriers=n, counts=counts, prior=training.prior)


def max_product(left: sp.csr_matrix, right: sp.csr_matrix, return_argmax: bool = False):
    """``out[i, j] = max over m of left[i, m] · right[m, j]`` (non-negative, sparse).

    The max-times product behind every projection here: with ``left`` a 0/1
    protein x unit incidence and ``right`` = ``θ`` it is "best unit wins";
    with ``left`` a score matrix and ``right`` a term -> ancestor-or-self
    indicator it is the max-propagation that makes scores hierarchy-consistent.
    With ``return_argmax`` also returns the arg-max ``m`` (as a CSR matrix of
    the same pattern; ``m`` stored as ``m + 1`` so zero stays "absent").
    """
    l_ptr, l_idx, l_val = left.indptr, left.indices, left.data
    r_ptr, r_idx, r_val = right.indptr, right.indices, right.data
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []
    args: list[np.ndarray] = []
    for i in range(left.shape[0]):
        mids = l_idx[l_ptr[i] : l_ptr[i + 1]]
        if mids.size == 0:
            continue
        weights = l_val[l_ptr[i] : l_ptr[i + 1]]
        starts = r_ptr[mids]
        lengths = r_ptr[mids + 1] - starts
        total = int(lengths.sum())
        if total == 0:
            continue
        # Concatenated ranges [start, start + length) without a Python loop.
        offsets = np.repeat(starts - np.cumsum(np.r_[0, lengths[:-1]]), lengths)
        gather = offsets + np.arange(total)
        c = r_idx[gather]
        v = r_val[gather] * np.repeat(weights, lengths)
        order = np.lexsort((-v, c))  # by column, best value first
        c, v = c[order], v[order]
        first = np.flatnonzero(np.r_[True, c[1:] != c[:-1]])
        rows.append(np.full(first.size, i))
        cols.append(c[first])
        vals.append(v[first])
        if return_argmax:
            args.append(np.repeat(mids, lengths)[order][first] + 1)
    shape = (left.shape[0], right.shape[1])
    if not rows:
        empty = sp.csr_matrix(shape)
        return (empty, empty.copy()) if return_argmax else empty
    r, c = np.concatenate(rows), np.concatenate(cols)
    out = sp.csr_matrix((np.concatenate(vals), (r, c)), shape=shape)
    if not return_argmax:
        return out
    return out, sp.csr_matrix(
        (np.concatenate(args).astype(np.float64), (r, c)), shape=shape
    )


def max_project(query_incidence: sp.csr_matrix, theta: sp.csr_matrix) -> sp.csr_matrix:
    """``score(q, t) = max over q's units u of θ(u, t)`` (sparse; absent = 0)."""
    return max_product(query_incidence, theta)


def ancestor_matrix(terms: Sequence[str], get_ancestors) -> sp.csr_matrix:
    """Term x term indicator: ``M[c, a] = 1`` when ``a`` is ``c`` or an ancestor of it.

    Only ancestors inside ``terms`` are kept (aspect roots never are).
    """
    index = {t: i for i, t in enumerate(terms)}
    rows, cols = [], []
    for i, t in enumerate(terms):
        rows.append(i)
        cols.append(i)
        for a in get_ancestors(t):
            j = index.get(a)
            if j is not None:
                rows.append(i)
                cols.append(j)
    n = len(terms)
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))


def make_consistent(scores: sp.csr_matrix, ancestors: sp.csr_matrix) -> sp.csr_matrix:
    """Raise every ancestor to the best score among its scored descendants."""
    return max_product(scores.tocsr(), ancestors)


def _fit_term_block(X, Y_block, Xq, C):
    from sklearn.linear_model import LogisticRegression

    out = np.zeros((Xq.shape[0], Y_block.shape[1]), dtype=np.float32)
    for j in range(Y_block.shape[1]):
        y = Y_block[:, j].toarray().ravel() > 0
        model = LogisticRegression(C=C, solver="liblinear", max_iter=200)
        model.fit(X, y)
        out[:, j] = model.predict_proba(Xq)[:, 1]
    return out


def lr_project(
    train_incidence: sp.csr_matrix,
    training: AspectTraining,
    query_incidence: sp.csr_matrix,
    min_positives: int = 5,
    min_carriers: int = 2,
    C: float = 1.0,
    n_jobs: int = 32,
    ancestors: sp.csr_matrix | None = None,
) -> sp.csr_matrix:
    """LR-InterPro (GOLabeler, Zhu et al. 2018): one logistic regression per term.

    Features are the binary unit incidences (units with ``≥ min_carriers``
    annotated carriers); terms need ``≥ min_positives`` annotated proteins.
    Unlike "best unit wins", a per-term model can weigh several units at once
    and learn from units that argue *against* a term. Output is
    max-propagated up the DAG when ``ancestors`` is given (per-term models
    are not hierarchy-consistent on their own).
    """
    from joblib import Parallel, delayed

    keep_units = np.flatnonzero(
        np.asarray(train_incidence.sum(axis=0)).ravel() >= min_carriers
    )
    X = train_incidence[:, keep_units].tocsr()
    Xq = query_incidence[:, keep_units].tocsr()
    pos = np.asarray(training.labels.sum(axis=0)).ravel()
    terms = np.flatnonzero(pos >= min_positives)
    Y = training.labels.tocsc()[:, terms]
    blocks = np.array_split(np.arange(terms.size), max(1, min(n_jobs * 4, terms.size)))
    parts = Parallel(n_jobs=n_jobs)(
        delayed(_fit_term_block)(X, Y[:, b], Xq, C) for b in blocks if b.size
    )
    dense = np.concatenate(parts, axis=1) if parts else np.zeros((Xq.shape[0], 0))
    # Proteins without a kept unit get only the intercept: no evidence, so 0.
    has_unit = np.diff(Xq.indptr) > 0
    dense[~has_unit] = 0.0
    rows, cols = np.nonzero(dense > 1e-4)
    out = sp.csr_matrix(
        (dense[rows, cols].astype(np.float64), (rows, terms[cols])),
        shape=(Xq.shape[0], len(training.terms)),
    )
    return make_consistent(out, ancestors) if ancestors is not None else out


# --------------------------------------------------------------------------- #
# Literature neighbours                                                        #
# --------------------------------------------------------------------------- #
def knn_project(
    query_emb: np.ndarray,
    train_emb: np.ndarray,
    train_labels: sp.csr_matrix,
    k: int = 25,
    temperature: float = 0.05,
    chunk: int = 2048,
    exclude_self: Sequence[int] | None = None,
    max_sim: float = 0.9999,
) -> sp.csr_matrix:
    """Similarity-weighted vote of each query's ``k`` nearest training proteins.

    Rows of both embedding matrices must be L2-normalised. Weights are
    ``exp((cos − 1) / temperature)``, so the vote leans on the closest
    neighbours; ``score(q, t)`` is the weighted fraction of neighbours with
    ``t``. ``exclude_self[i]`` (a training row index, or -1) is removed from
    query ``i``'s neighbours — for cross-validation with shared rows.
    Neighbours with ``cos ≥ max_sim`` are dropped too: an identical vector is
    the same gene under another accession, an identifier artefact rather than
    knowledge to project.
    """
    n_q = query_emb.shape[0]
    k = min(k, train_emb.shape[0])
    out_rows: list[sp.csr_matrix] = []
    for start in range(0, n_q, chunk):
        q = query_emb[start : start + chunk]
        sim = q @ train_emb.T
        if exclude_self is not None:
            ex = np.asarray(exclude_self[start : start + chunk])
            ok = ex >= 0
            sim[np.flatnonzero(ok), ex[ok]] = -np.inf
        sim[sim >= max_sim] = -np.inf
        nn = np.argpartition(-sim, k - 1, axis=1)[:, :k]
        nn_sim = np.take_along_axis(sim, nn, axis=1)
        w = np.exp((nn_sim - 1.0) / temperature)
        w[~np.isfinite(nn_sim)] = 0.0
        w = w.astype(np.float64)
        w /= np.maximum(w.sum(axis=1, keepdims=True), 1e-300)
        weights = sp.csr_matrix(
            (w.ravel(), (np.repeat(np.arange(q.shape[0]), k), nn.ravel())),
            shape=(q.shape[0], train_emb.shape[0]),
        )
        out_rows.append((weights @ train_labels).tocsr())
    return (
        sp.vstack(out_rows).tocsr()
        if out_rows
        else sp.csr_matrix((0, train_labels.shape[1]))
    )


def load_literature(
    path: Path, min_papers: int = 1, center: bool = True
) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    """``({accession: unit vector}, {accession: n papers})`` from
    ``scripts/build_literature_features.py`` output."""
    data = np.load(path, allow_pickle=False)
    accs = data["accessions"].astype(str)
    emb = data["emb"].astype(np.float32)
    n = data["n_papers"].astype(int)
    keep = n >= min_papers
    if center:
        # Averaging hundreds of papers pulls every gene towards one shared
        # "biomedical" direction (random genes sit at cos ≈ 0.9). Removing the
        # mean (no labels involved) restores the spread neighbours need.
        emb = emb - emb[keep].mean(axis=0)
        emb /= np.maximum(np.linalg.norm(emb, axis=1, keepdims=True), 1e-12)
    return (
        {a: e for a, e, k in zip(accs, emb, keep) if k},
        {a: int(c) for a, c, k in zip(accs, n, keep) if k},
    )


def load_h5_embeddings(
    paths: Sequence[Path], accessions: set[str] | None = None, center: bool = True
) -> dict[str, np.ndarray]:
    """``{accession: unit vector}`` from UniProt's per-protein ProtT5 HDF5 files.

    A sequence-only representation (no function labels went into ProtT5), so a
    neighbour vote in this space projects knowledge to any protein with a
    sequence — domains or not. The first file listing an accession wins.
    Centred and re-normalised like :func:`load_literature`.
    """
    import h5py

    out: dict[str, np.ndarray] = {}
    for path in paths:
        with h5py.File(path, "r") as handle:
            for acc in handle.keys():
                if acc in out or (accessions is not None and acc not in accessions):
                    continue
                out[acc] = np.asarray(handle[acc], dtype=np.float32)
    if not out:
        return out
    keys = list(out)
    emb = np.stack([out[k] for k in keys])
    emb /= np.maximum(np.linalg.norm(emb, axis=1, keepdims=True), 1e-12)
    if center:
        emb -= emb.mean(axis=0)
        emb /= np.maximum(np.linalg.norm(emb, axis=1, keepdims=True), 1e-12)
    return dict(zip(keys, emb))


def load_literature_embeddings(
    path: Path, min_papers: int = 1
) -> dict[str, np.ndarray]:
    """``{accession: unit vector}`` (see :func:`load_literature`)."""
    return load_literature(path, min_papers)[0]


# --------------------------------------------------------------------------- #
# Output                                                                       #
# --------------------------------------------------------------------------- #
def to_prediction_dict(
    scores: sp.csr_matrix,
    proteins: Sequence[str],
    terms: Sequence[str],
    min_score: float = 0.0,
) -> dict[str, dict[str, float]]:
    """``{protein: {term: score}}`` from a query x term score matrix."""
    scores = scores.tocsr()
    out: dict[str, dict[str, float]] = {}
    for i, p in enumerate(proteins):
        lo, hi = scores.indptr[i], scores.indptr[i + 1]
        if hi == lo:
            continue
        idx, val = scores.indices[lo:hi], scores.data[lo:hi]
        keep = val > min_score
        if keep.any():
            out[p] = {terms[j]: float(v) for j, v in zip(idx[keep], val[keep])}
    return out


def weighted_mean(
    mats: Sequence[sp.csr_matrix],
    weights: Sequence[float],
    skip_missing: bool = False,
) -> sp.csr_matrix:
    """Weighted mean of channel score matrices (one shape, one term order).

    By default a channel with no score for a protein counts 0, so a protein
    missing a channel is capped at the other channels' share of the weight
    (the rule frozen on dev, ``KNOWLEDGE_PROJECTION.md`` §1.6). With
    ``skip_missing`` a channel with no score at all for a protein (no
    papers, no embedding, no unit with annotated carriers) leaves that
    protein's denominator instead; on dev that raised pooled AUPRC slightly
    and lowered the full arm's KY50, so it is a sensitivity option only.
    Hierarchy consistency is kept either way: each channel is consistent and
    a protein's weights are the same for all its terms.
    """
    if len(mats) != len(weights) or not mats:
        raise ValueError("one weight per channel")
    n = mats[0].shape[0]
    total = sp.csr_matrix(mats[0].shape)
    denom = np.zeros(n)
    for m, w in zip(mats, weights):
        m = m.tocsr()
        total = total + w * m
        denom += w * (np.diff(m.indptr) > 0) if skip_missing else w
    scale = np.where(denom > 0, 1.0 / np.maximum(denom, 1e-300), 0.0)
    return (sp.diags(scale) @ total).tocsr()


def noisy_or(a: sp.csr_matrix, b: sp.csr_matrix) -> sp.csr_matrix:
    """``1 − (1 − a)(1 − b)`` elementwise on two score matrices of one shape."""
    return (a + b - a.multiply(b)).tocsr()
