"""End-to-end test of the surprise-score driver on a synthetic world.

Builds a tiny protein universe with three planted situations and checks the
ranking separates them:

* an **emergent** pair — neither domain predicts the term alone, the combination
  always does;
* a **redundant-signature** pair — two signatures over the same region, so the
  "combination" is one domain described twice;
* an **explained** pair — one constituent already predicts the term.
"""

import importlib.util
from pathlib import Path

import pytest

from src.domain_annotation_parser import DomainAnnotation, ProteinDomainArchitecture
from src.surprise_score import (
    EmergenceEvidence,
    apply_fdr,
    conditional_rate,
    proper_subfeatures,
    score_candidate,
)

DRIVER_PATH = (
    Path(__file__).resolve().parents[2] / "scripts" / "rank_surprising_associations.py"
)


@pytest.fixture(scope="module")
def driver():
    """Import the CLI driver as a module so its helpers can be exercised."""
    spec = importlib.util.spec_from_file_location("rank_surprising", DRIVER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Each scenario gets its own domains so their rates stay independent, and the
# coordinates fix the order in which a supra-domain id is spelled (features are
# built from the positionally sorted domain list).
#   A + B    — separate regions, the emergent pair
#   D + D2   — two signatures over one region, the redundant pair
#   C + E    — C already predicts the term on its own
#   F        — filler, so the term's background rate stays low
REGIONS = {
    "A": (10, 100),
    "B": (200, 300),
    "C": (400, 500),
    "D": (600, 700),
    "D2": (605, 695),
    "E": (800, 900),
    "F": (1000, 1100),
}


def architecture(protein_id, domain_ids):
    annotations = [
        DomainAnnotation(
            protein_id=protein_id,
            interpro_id=domain_id,
            interpro_name=domain_id,
            signature_id=f"SIG{domain_id}",
            start=REGIONS[domain_id][0],
            end=REGIONS[domain_id][1],
        )
        for domain_id in domain_ids
    ]
    annotations.sort(key=lambda a: a.start)
    ordered = [a.interpro_id for a in annotations]
    supra = [
        ",".join(ordered[i : i + length])
        for length in (2, 3)
        for i in range(len(ordered) - length + 1)
    ]
    return ProteinDomainArchitecture(
        protein_id=protein_id,
        single_domains=ordered,
        supra_domains=supra,
        domain_annotations=annotations,
    )


@pytest.fixture
def world():
    """A 190-protein universe with the three planted situations."""
    architectures = {}
    annotations = {}

    def add(prefix, count, domains, annotated):
        for i in range(count):
            pid = f"{prefix}{i}"
            architectures[pid] = architecture(pid, domains)
            annotations[pid] = {"T"} if annotated else {"OTHER"}

    add("EMERGENT", 4, ["A", "B"], True)  # only the pair carries the term
    add("ONLYA", 60, ["A"], False)  # A alone: silent
    add("ONLYB", 60, ["B"], False)  # B alone: silent
    add("REDUNDANT", 4, ["D", "D2"], True)  # one region, two signatures
    add("ONLYD", 30, ["D"], False)
    add("EXPLAINED", 4, ["C", "E"], True)  # C already predicts the term
    add("ONLYC", 8, ["C"], True)
    add("ONLYE", 20, ["E"], False)
    return architectures, annotations


def score_all(driver, world, candidates):
    """Run the driver's indexing + scoring path over the synthetic world."""
    architectures, annotations = world
    universe = set(architectures)
    wanted = set()
    subfeatures = {}
    for feature, parts in candidates:
        subfeatures[feature] = proper_subfeatures(parts)
        wanted.add(feature)
        wanted.update(parts)
        wanted.update(subfeatures[feature])
    index = driver.build_feature_index(architectures, universe, wanted)

    annotated = {p for p, terms in annotations.items() if "T" in terms}
    background = len(annotated) / len(universe)

    results = []
    for feature, parts in candidates:
        carriers = index.get(feature, set())
        rates = [
            conditional_rate(
                len(index[part] & annotated), len(index[part]), background, 1.0
            )
            for part in dict.fromkeys(parts)
            if index.get(part)
        ]
        evidence = EmergenceEvidence(
            feature=feature,
            term="T",
            n_feature=len(carriers),
            n_both=len(carriers & annotated),
            single_rates=tuple(rates),
            part_rates=(),
            background_rate=background,
            q_value=1e-6,
        )
        overlap = driver.measure_region_overlap(
            feature,
            parts,
            carriers & annotated,
            architectures,
            # This world's architectures are keyed by InterPro entry, which is
            # what DomainAnnotationParser(domain_key="interpro") reports.
            lambda annotation: annotation.interpro_id,
        )
        results.append(score_candidate(evidence, overlap))
    return {r.feature: r for r in apply_fdr(results)}


@pytest.fixture
def scored(driver, world):
    return score_all(
        driver,
        world,
        [("A,B", ("A", "B")), ("D,D2", ("D", "D2")), ("C,E", ("C", "E"))],
    )


class TestSurpriseRanking:
    def test_emergent_pair_beats_its_parts(self, scored):
        emergent = scored["A,B"]
        assert emergent.n_both == 4
        assert emergent.observed_rate == 1.0
        assert emergent.expected_rate < 0.3
        assert emergent.q_emergence < 0.05

    def test_redundant_signatures_are_detected_by_region_overlap(self, scored):
        assert scored["D,D2"].region_overlap > 0.9
        assert scored["D,D2"].distinctness < 0.1

    def test_redundant_pair_scores_below_the_emergent_one(self, scored):
        assert scored["A,B"].surprise > scored["D,D2"].surprise

    def test_pair_explained_by_one_constituent_is_not_surprising(self, scored):
        explained = scored["C,E"]
        # C alone predicts T, so the expectation is already high.
        assert explained.expected_rate > 0.5
        assert explained.surprise < scored["A,B"].surprise

    def test_uninformative_constituents_are_counted(self, scored):
        # A and B each predict the term below its background rate; C does not.
        assert scored["A,B"].uninformative_constituents == 2
        assert scored["C,E"].uninformative_constituents == 0


class TestAssociationsFile:
    def test_reads_supra_rows_and_skips_singles(self, driver, tmp_path):
        path = tmp_path / "assoc.tsv"
        path.write_text(
            "domain\tgo_term\tp_value\tadj_p_value\todds_ratio\thyper_score\t"
            "domain_type\tconstituent_domains\tn_observations\n"
            "IPR1\tGO:1\t1e-9\t1e-6\tinf\t100\tsingle\t-\t10\n"
            "IPR1,IPR2\tGO:1\t1e-9\t2e-6\tinf\t100\tsupra_pair\tIPR1,IPR2\t4\n"
        )
        rows = driver.parse_associations(path, "go_term")
        assert rows == [("IPR1,IPR2", "GO:1", 2e-6, ("IPR1", "IPR2"))]

    def test_wrong_term_column_fails_loudly(self, driver, tmp_path):
        path = tmp_path / "assoc.tsv"
        path.write_text("domain\tgo_term\tadj_p_value\tconstituent_domains\n")
        with pytest.raises(SystemExit, match="ec_term"):
            driver.parse_associations(path, "ec_term")
