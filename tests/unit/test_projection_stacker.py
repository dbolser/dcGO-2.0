"""Tests for the stacked projection features (src/projection_stacker.py)."""

import numpy as np
import pytest

from src import knowledge_projection as kp
from src import projection_stacker as ps

PARENTS = {
    "GO:a": {"GO:b", "GO:0008150"},
    "GO:b": {"GO:0008150"},
    "GO:c": {"GO:0008150"},
}
ASPECT = {t: "BP" for t in ("GO:a", "GO:b", "GO:c", "GO:0008150")}


def anc(t):
    return PARENTS.get(t, set())


def _training():
    ann = {f"P{i}": ({"GO:a"} if i % 2 else {"GO:c"}) for i in range(8)}
    return kp.build_aspect_training(ann, anc, ASPECT, "BP")


def test_features_align_with_theta_and_support():
    tr = _training()
    fam = {
        "ipr": {**{f"P{i}": {"D1" if i % 2 else "D2"} for i in range(8)}, "Q": {"D1"}},
        "supra": {},
        "arch": {},
        "sig": {},
    }
    feats = ps.build_features(tr, fam, ["Q"], strength=1.0, ic={"GO:a": 1.0})
    names = list(ps.FEATURES)
    by_term = {tr.terms[c]: feats.X[i] for i, c in enumerate(feats.cols)}
    assert set(by_term) == {"GO:a", "GO:b"}  # D1 carriers only have a (and b)
    row = by_term["GO:a"]
    # D1: n = 4 carriers, k = 4 with GO:a, prior 0.5, s = 1 -> (4.5) / 5
    assert row[names.index("theta_max")] == pytest.approx(0.9)
    assert row[names.index("theta_ipr")] == pytest.approx(0.9)
    assert row[names.index("support_log_n")] == pytest.approx(np.log1p(4))
    assert row[names.index("support_k")] == pytest.approx(np.log1p(4), rel=1e-5)
    assert row[names.index("ic")] == pytest.approx(1.0)
    assert feats.labels({"Q": {"GO:a"}}).tolist() == [
        tr.terms[c] == "GO:a" for c in feats.cols
    ]


def test_consistency_lifts_parents():
    terms = ("GO:a", "GO:b")
    m = kp.ancestor_matrix(terms, anc)
    import scipy.sparse as sp

    s = sp.csr_matrix(np.array([[0.8, 0.3]]))
    np.testing.assert_allclose(kp.make_consistent(s, m).toarray(), [[0.8, 0.8]])


def test_cross_fit_is_out_of_fold():
    pytest.importorskip("sklearn")
    rng = np.random.default_rng(0)
    n = 4000
    X = rng.random((n, len(ps.FEATURES))).astype(np.float32)
    y = X[:, 0] + 0.1 * rng.standard_normal(n) > 0.6
    feats = ps.CandidateFeatures(
        tuple(f"q{i}" for i in range(400)),
        ("t",),
        np.arange(n) % 400,
        np.zeros(n, int),
        X,
    )
    p = ps.cross_fit(feats, y, n_folds=2)
    assert p.shape == (n,) and (p > 0).all()
    assert np.corrcoef(p, X[:, 0])[0, 1] > 0.8
