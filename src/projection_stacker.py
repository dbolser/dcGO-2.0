"""Learned combination of projection evidence (GOLabeler-style stacking).

:mod:`src.knowledge_projection` gives each (protein, term) several pieces of
evidence — the best unit's ``θ``, which family that unit came from, how many
annotated proteins stand behind it, how many of the protein's units agree,
the literature neighbours' vote. A fixed rule (max, noisy-OR) has to guess
how they trade off; this module learns it, the way GOLabeler (Zhu et al.
2018) and NetGO learn-to-rank over component predictors.

The model is a gradient-boosted classifier on one row per candidate
(protein, term). It is **monotone increasing** in every evidence score, so
more evidence can never lower a prediction; that keeps it from learning
quirks of a small training cohort. Its output is made hierarchy-consistent
by max-propagation, so a parent never scores below a child.

Training rows must come from a split that ends before the one predicted —
the dev split's no-knowledge cohort (t0 features → t1 truth) for the final
split, never the final split itself.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from src import knowledge_projection as kp

#: Unit families with their own max-θ feature.
FAMILIES: tuple[str, ...] = ("ipr", "supra", "arch", "sig")

#: Feature name -> whether the model must be monotone increasing in it.
MONOTONE: dict[str, int] = {
    "theta_max": 1,
    "theta_ipr": 1,
    "theta_supra": 1,
    "theta_arch": 1,
    "theta_sig": 1,
    "lift_max": 1,
    "support_log_n": 0,
    "support_k": 1,
    "n_agree": 1,
    "frac_agree": 1,
    "lit_sharp": 1,
    "lit_smooth": 1,
    "lit_top1_sim": 0,
    "log_n_papers": 0,
    "log_n_units": 0,
    "log_prior": 0,
    "ic": 0,
}
FEATURES: tuple[str, ...] = tuple(MONOTONE)


@dataclass(frozen=True)
class CandidateFeatures:
    """One row per candidate (query, term) for one aspect."""

    queries: tuple[str, ...]
    terms: tuple[str, ...]
    rows: np.ndarray  # query index
    cols: np.ndarray  # term index into ``terms``
    X: np.ndarray  # float32, len(rows) x len(FEATURES)

    def labels(self, truth: Mapping[str, set[str]]) -> np.ndarray:
        """1 where the candidate term is in its query's truth set."""
        return np.fromiter(
            (
                self.terms[c] in truth.get(self.queries[r], ())
                for r, c in zip(self.rows, self.cols)
            ),
            dtype=bool,
            count=self.rows.size,
        )

    def subset(self, keep: np.ndarray) -> CandidateFeatures:
        return CandidateFeatures(
            self.queries, self.terms, self.rows[keep], self.cols[keep], self.X[keep]
        )

    def to_matrix(self, scores: np.ndarray) -> sp.csr_matrix:
        return sp.csr_matrix(
            (scores.astype(np.float64), (self.rows, self.cols)),
            shape=(len(self.queries), len(self.terms)),
        )


def _values_at(mat: sp.csr_matrix, keys: np.ndarray, n_cols: int) -> np.ndarray:
    """``mat`` read at flattened ``row * n_cols + col`` keys (0 where absent)."""
    coo = mat.tocoo()
    mkeys = coo.row.astype(np.int64) * n_cols + coo.col
    order = np.argsort(mkeys)
    mkeys, mvals = mkeys[order], coo.data[order]
    pos = np.searchsorted(mkeys, keys)
    pos_c = np.minimum(pos, max(mkeys.size - 1, 0))
    hit = (
        (pos < mkeys.size) & (mkeys[pos_c] == keys)
        if mkeys.size
        else np.zeros(keys.size, bool)
    )
    out = np.zeros(keys.size)
    out[hit] = mvals[pos_c[hit]]
    return out


def build_features(
    training: kp.AspectTraining,
    families: Mapping[str, Mapping[str, set[str]]],
    queries: Sequence[str],
    strength: float,
    ic: Mapping[str, float],
    embeddings: Mapping[str, np.ndarray] | None = None,
    n_papers: Mapping[str, int] | None = None,
    knn_k: int = 25,
    temperatures: tuple[float, float] = (0.05, 0.2),
) -> CandidateFeatures:
    """Every evidence feature for ``queries`` against ``training``."""
    queries = tuple(queries)
    n_t = len(training.terms)
    merged = kp.merge_units(*(families[f] for f in FAMILIES if f in families))
    vocab = kp.UnitVocabulary(merged)
    evidence = kp.unit_term_evidence(vocab.incidence(training.proteins), training)
    theta = evidence.posterior(strength)
    q_inc = vocab.incidence(queries)

    theta_max, arg = kp.max_product(q_inc, theta, return_argmax=True)
    per_family = {}
    for fam in FAMILIES:
        units = {u for us in families.get(fam, {}).values() for u in us}
        mask = np.zeros(len(vocab.units))
        mask[[vocab.index[u] for u in units if u in vocab.index]] = 1.0
        fam_inc = (q_inc @ sp.diags(mask)).tocsr()
        fam_inc.eliminate_zeros()
        per_family[fam] = kp.max_product(fam_inc, theta)
    lift = theta.multiply(1.0 / np.maximum(training.prior, 1e-9)).tocsr()
    lift_max = kp.max_product(q_inc, lift)
    agree = (q_inc @ (evidence.counts > 0).astype(np.float64)).tocsr()

    lit = []
    top1 = np.zeros(len(queries))
    if embeddings:
        lit, top1 = _literature(training, embeddings, queries, knn_k, temperatures)

    pattern = (theta_max > 0).astype(np.float64)
    for m in lit:
        pattern = pattern + (m > 0)
    pattern = pattern.tocoo()
    rows, cols = pattern.row.astype(np.int64), pattern.col.astype(np.int64)
    keys = rows * n_t + cols

    # The argmax unit's support: n carriers, and k of them with the term.
    arg_unit = _values_at(arg, keys, n_t).astype(np.int64) - 1
    has_unit = arg_unit >= 0
    n_sup = np.where(has_unit, evidence.n_carriers[np.maximum(arg_unit, 0)], 0.0)
    th = _values_at(theta_max, keys, n_t)
    k_sup = np.where(
        has_unit, th * (n_sup + strength) - strength * training.prior[cols], 0
    )

    n_units = np.diff(q_inc.indptr)[rows]
    n_agree = _values_at(agree, keys, n_t)
    term_ic = np.array([ic.get(t, 0.0) for t in training.terms])[cols]
    q_papers = np.array([(n_papers or {}).get(q, 0) for q in queries])[rows]

    cols_out = {
        "theta_max": th,
        "theta_ipr": _values_at(per_family["ipr"], keys, n_t),
        "theta_supra": _values_at(per_family["supra"], keys, n_t),
        "theta_arch": _values_at(per_family["arch"], keys, n_t),
        "theta_sig": _values_at(per_family["sig"], keys, n_t),
        "lift_max": np.log1p(_values_at(lift_max, keys, n_t)),
        "support_log_n": np.log1p(n_sup),
        "support_k": np.log1p(np.maximum(k_sup, 0)),
        "n_agree": np.log1p(n_agree),
        "frac_agree": n_agree / np.maximum(n_units, 1),
        "lit_sharp": _values_at(lit[0], keys, n_t) if lit else np.zeros(keys.size),
        "lit_smooth": _values_at(lit[1], keys, n_t) if lit else np.zeros(keys.size),
        "lit_top1_sim": top1[rows],
        "log_n_papers": np.log1p(q_papers),
        "log_n_units": np.log1p(n_units),
        "log_prior": np.log10(np.maximum(training.prior[cols], 1e-6)),
        "ic": term_ic,
    }
    X = np.column_stack([cols_out[f] for f in FEATURES]).astype(np.float32)
    return CandidateFeatures(queries, training.terms, rows, cols, X)


def _literature(training, embeddings, queries, k, temperatures):
    train_rows = [i for i, p in enumerate(training.proteins) if p in embeddings]
    q_rows = [i for i, q in enumerate(queries) if q in embeddings]
    n_q, n_t = len(queries), len(training.terms)
    if not train_rows or not q_rows:
        empty = sp.csr_matrix((n_q, n_t))
        return [empty, empty.copy()], np.zeros(n_q)
    train_emb = np.stack([embeddings[training.proteins[i]] for i in train_rows])
    q_emb = np.stack([embeddings[queries[i]] for i in q_rows])
    labels = training.labels[train_rows]
    expand = sp.csr_matrix(
        (np.ones(len(q_rows)), (q_rows, np.arange(len(q_rows)))),
        shape=(n_q, len(q_rows)),
    )
    mats = [
        (expand @ kp.knn_project(q_emb, train_emb, labels, k=k, temperature=t)).tocsr()
        for t in temperatures
    ]
    top1 = np.zeros(n_q)
    for start in range(0, len(q_rows), 4096):
        sim = q_emb[start : start + 4096] @ train_emb.T
        sim[sim >= 0.9999] = -np.inf  # same gene, see knn_project
        top1[q_rows[start : start + 4096]] = sim.max(axis=1)
    return mats, top1


# --------------------------------------------------------------------------- #
# Model                                                                        #
# --------------------------------------------------------------------------- #
def fit(
    X: np.ndarray,
    y: np.ndarray,
    use: Sequence[str] = FEATURES,
    seed: int = 0,
    max_iter: int = 300,
):
    """Monotone gradient-boosted classifier on candidate rows (all of FEATURES)."""
    from sklearn.ensemble import HistGradientBoostingClassifier

    cols = [FEATURES.index(f) for f in use]
    model = HistGradientBoostingClassifier(
        learning_rate=0.05,
        max_iter=max_iter,
        max_leaf_nodes=31,
        min_samples_leaf=200,
        l2_regularization=1.0,
        monotonic_cst=[MONOTONE[f] for f in use],
        early_stopping=True,
        validation_fraction=0.1,
        n_iter_no_change=20,
        random_state=seed,
    )
    model.fit(X[:, cols], y)
    model.feature_names_used_ = tuple(use)
    return model


def cross_fit(
    features: CandidateFeatures,
    y: np.ndarray,
    use: Sequence[str] = FEATURES,
    n_folds: int = 2,
    seed: int = 0,
) -> np.ndarray:
    """Out-of-fold probabilities, folds drawn over queries (never over rows)."""
    rng = np.random.default_rng(seed)
    fold_of_query = rng.integers(0, n_folds, size=len(features.queries))
    fold = fold_of_query[features.rows]
    out = np.zeros(features.rows.size)
    cols = [FEATURES.index(f) for f in use]
    for f in range(n_folds):
        train, test = fold != f, fold == f
        if not test.any() or not y[train].any():
            continue
        model = fit(features.X[train], y[train], use, seed=seed)
        out[test] = model.predict_proba(features.X[test][:, cols])[:, 1]
    return out


def predict(
    model, features: CandidateFeatures, ancestors: sp.csr_matrix | None = None
) -> sp.csr_matrix:
    """Query x term probabilities, max-propagated when ``ancestors`` is given."""
    cols = [FEATURES.index(f) for f in model.feature_names_used_]
    p = model.predict_proba(features.X[:, cols])[:, 1]
    return finish(features, p, ancestors)


def finish(
    features: CandidateFeatures, p: np.ndarray, ancestors: sp.csr_matrix | None
) -> sp.csr_matrix:
    """Probabilities as a query x term matrix, max-propagated up the DAG."""
    scores = features.to_matrix(np.maximum(p, 1e-12))
    return kp.make_consistent(scores, ancestors) if ancestors is not None else scores
