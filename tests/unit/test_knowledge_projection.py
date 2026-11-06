"""Tests for calibrated unit projection (src/knowledge_projection.py)."""

import numpy as np
import pytest
import scipy.sparse as sp

from src import knowledge_projection as kp

# Tiny DAG: leaf -> mid -> BP root; other is a separate BP child of the root.
PARENTS = {
    "GO:leaf": {"GO:mid", "GO:0008150"},
    "GO:mid": {"GO:0008150"},
    "GO:other": {"GO:0008150"},
    "GO:mf": {"GO:0003674"},
}
ASPECT = {
    "GO:leaf": "BP",
    "GO:mid": "BP",
    "GO:other": "BP",
    "GO:0008150": "BP",
    "GO:mf": "MF",
    "GO:0003674": "MF",
}


def ancestors(t):
    return PARENTS.get(t, set())


def test_training_is_aspect_annotated_and_propagated():
    tr = kp.build_aspect_training(
        {"P1": {"GO:leaf"}, "P2": {"GO:mf"}, "P3": {"GO:other"}},
        ancestors,
        ASPECT,
        "BP",
    )
    # P2 has no BP annotation, so it says nothing about BP terms.
    assert tr.proteins == ("P1", "P3")
    assert set(tr.terms) == {"GO:leaf", "GO:mid", "GO:other"}  # root dropped
    i = tr.term_index
    assert tr.prior[i["GO:mid"]] == pytest.approx(0.5)


def test_posterior_shrinks_to_prior_and_max_projects():
    tr = kp.build_aspect_training(
        {f"P{i}": ({"GO:leaf"} if i < 3 else {"GO:other"}) for i in range(4)},
        ancestors,
        ASPECT,
        "BP",
    )
    units = {"P0": {"D1"}, "P1": {"D1"}, "P2": {"D1", "D2"}, "P3": {"D2"}, "Q": {"D1"}}
    vocab = kp.UnitVocabulary(units)
    ev = kp.unit_term_evidence(vocab.incidence(tr.proteins), tr)
    theta = ev.posterior(strength=1.0)
    i, d1 = tr.term_index, vocab.index["D1"]
    prior_leaf = 3 / 4
    assert theta[d1, i["GO:leaf"]] == pytest.approx((3 + prior_leaf) / (3 + 1))
    q = kp.max_project(vocab.incidence(["Q"]), theta)
    assert q[0, i["GO:leaf"]] == pytest.approx(theta[d1, i["GO:leaf"]])
    # Parent never below child.
    assert q[0, i["GO:mid"]] >= q[0, i["GO:leaf"]]
    assert q[0, i["GO:other"]] == 0.0  # D1 never seen with it


def test_max_project_takes_best_unit():
    theta = sp.csr_matrix(np.array([[0.2, 0.0], [0.6, 0.1]]))
    inc = sp.csr_matrix(np.array([[1, 1], [1, 0], [0, 0]], dtype=float))
    out = kp.max_project(inc, theta).toarray()
    np.testing.assert_allclose(out, [[0.6, 0.1], [0.2, 0.0], [0.0, 0.0]])


def test_knn_votes_with_nearest_and_excludes_self():
    train = np.eye(3, dtype=np.float32)
    labels = sp.csr_matrix(np.array([[1, 0], [0, 1], [0, 1]], dtype=float))
    q = np.array([[1, 0.1, 0]], dtype=np.float32)
    q /= np.linalg.norm(q)
    out = kp.knn_project(q, train, labels, k=1).toarray()
    np.testing.assert_allclose(out, [[1.0, 0.0]])
    out = kp.knn_project(q, train, labels, k=1, exclude_self=[0]).toarray()
    assert out[0, 1] == pytest.approx(1.0)


def test_knn_drops_identical_vectors():
    """Same gene under another accession: not a neighbour."""
    train = np.eye(3, dtype=np.float32)
    labels = sp.csr_matrix(np.array([[1, 0], [0, 1], [0, 1]], dtype=float))
    q = np.array([[1, 0, 0]], dtype=np.float32)
    out = kp.knn_project(q, train, labels, k=1).toarray()
    assert out[0, 0] == 0.0


def test_units_helpers():
    doms = {"P": ["A", "B", "C", "D", "A,B"]}
    assert kp.architecture_units(doms) == {"P": {"ARCH:A+B+C+D"}}
    supra = kp.bag_supra_units({"P": ["B", "A", "C"]})["P"]
    assert supra == {"A,B", "A,C", "B,C", "A,B,C"}


def test_noisy_or():
    a = sp.csr_matrix(np.array([[0.5, 0.0]]))
    b = sp.csr_matrix(np.array([[0.5, 0.2]]))
    np.testing.assert_allclose(kp.noisy_or(a, b).toarray(), [[0.75, 0.2]])


def test_read_restricted_units(tmp_path):
    import gzip

    rows = [
        ("X1", "IPR1", "n", "PF1", 1, 50),
        ("X1", "IPR2", "n", "PTHR9:SF2", 60, 120),
        ("X1", "IPR3", "n", "SM3", 130, 135),  # too short for an entry
        ("X1", "IPR4", "n", "PF4", 140, 200),
        ("X1", "IPR5", "n", "PF5", 210, 300),
        ("X2", "IPR1", "n", "PF1", 1, 50),
        ("Y", "IPR1", "n", "PF1", 1, 50),
    ]
    path = tmp_path / "p2i.dat.gz"
    with gzip.open(path, "wt") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    vocab = {"IPR1", "IPR1,IPR2", "SIG:PTHR9", "ARCH:IPR1+IPR2+IPR4+IPR5", "IPR9"}
    out = kp.read_restricted_units(path, {"X1", "X2"}, vocab)
    assert out["X1"] == {"IPR1", "IPR1,IPR2", "SIG:PTHR9", "ARCH:IPR1+IPR2+IPR4+IPR5"}
    assert out["X2"] == {"IPR1"}
    assert "Y" not in out
    fams = kp.split_families(out)
    assert fams["supra"]["X1"] == {"IPR1,IPR2"} and fams["arch"]["X1"]


def test_lr_project_learns_a_unit_and_zeroes_unitless_queries():
    pytest.importorskip("sklearn")
    rng = np.random.default_rng(0)
    n = 60
    has_d1 = rng.random(n) < 0.5
    ann = {f"P{i}": ({"GO:leaf"} if has_d1[i] else {"GO:other"}) for i in range(n)}
    tr = kp.build_aspect_training(ann, ancestors, ASPECT, "BP")
    units = {f"P{i}": ({"D1"} if has_d1[i] else {"D2"}) for i in range(n)}
    units["Q1"] = {"D1"}
    units["Q2"] = {"D2"}
    vocab = kp.UnitVocabulary(units)
    anc = kp.ancestor_matrix(tr.terms, ancestors)
    out = kp.lr_project(
        vocab.incidence(tr.proteins),
        tr,
        vocab.incidence(["Q1", "Q2", "Q3"]),
        min_positives=3,
        n_jobs=1,
        ancestors=anc,
    ).toarray()
    i = tr.term_index
    assert out[0, i["GO:leaf"]] > 0.5 > out[1, i["GO:leaf"]]
    assert out[0, i["GO:mid"]] >= out[0, i["GO:leaf"]]  # consistent
    assert out[2].sum() == 0  # no units, no prediction


def test_vectorised_trainings_match_loop():
    ann = {
        "P1": {"GO:leaf"},
        "P2": {"GO:mf", "GO:other"},
        "P3": {"GO:other"},
        "P4": {"GO:mf"},
    }
    many = kp.build_aspect_trainings(
        ann, ancestors, ASPECT, ["BP", "MF"], exclude={"BP": ["P3"]}
    )
    for aspect, ex in (("BP", ["P3"]), ("MF", [])):
        one = kp.build_aspect_training(ann, ancestors, ASPECT, aspect, exclude=ex)
        got = many[aspect]
        assert got.proteins == one.proteins and got.terms == one.terms
        np.testing.assert_array_equal(got.labels.toarray(), one.labels.toarray())
        np.testing.assert_allclose(got.prior, one.prior)


def test_weighted_mean_skips_missing_channels():
    a = sp.csr_matrix(np.array([[0.8, 0.2], [0.6, 0.0]]))
    b = sp.csr_matrix(np.array([[0.4, 0.0], [0.0, 0.0]]))  # no data for row 1
    out = kp.weighted_mean([a, b], [1.0, 1.0], skip_missing=True).toarray()
    np.testing.assert_allclose(out, [[0.6, 0.1], [0.6, 0.0]])
    plain = kp.weighted_mean([a, b], [1.0, 1.0]).toarray()
    np.testing.assert_allclose(plain, [[0.6, 0.1], [0.3, 0.0]])


def test_knn_all_neighbours_excluded_gives_no_nan():
    train = np.array([[1.0, 0.0]], dtype=np.float32)
    labels = sp.csr_matrix(np.array([[1.0, 1.0]]))
    out = kp.knn_project(np.array([[1.0, 0.0]], dtype=np.float32), train, labels, k=5)
    assert not np.isnan(out.toarray()).any()
