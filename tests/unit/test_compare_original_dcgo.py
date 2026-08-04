"""Unit tests for the published-dcGO reference parsers (VALIDATION_PLAN §3).

These parsers are load-bearing: every agreement number in §3.1 is a join over
what they return, and the MySQL dump is read with a regex rather than a real SQL
parser, so its shape is worth pinning down.
"""

import gzip
import importlib.util
from pathlib import Path

import pytest

_MODULE_PATH = (
    Path(__file__).resolve().parents[2] / "validation" / "compare_original_dcgo.py"
)
_spec = importlib.util.spec_from_file_location("compare_original_dcgo", _MODULE_PATH)
compare = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(compare)


@pytest.fixture
def dump(tmp_path):
    """A miniature Domain2GO.sql.gz with the shapes the real dump contains."""
    path = tmp_path / "Domain2GO.sql.gz"
    sql = (
        "CREATE TABLE `GO_mapping` (\n  `id` mediumint(8) unsigned NOT NULL\n);\n"
        "INSERT INTO `GO_mapping` VALUES "
        # sf, significant, not published (empty inherited_from_all), no h-score
        "(52540,'sf',0003674,0.5,0.00000001,'','',NULL,NULL),"
        # sf, published (non-empty inherited_from_all, commas inside the quotes)
        "(52540,'sf',0005524,1,0.002,'','0042995,0009288',0,95.45),"
        # fa rows must not leak into the sf result
        "(54870,'fa',0009055,0.1,0.0001,'','',NULL,NULL),"
        # a GO id whose zerofill has leading zeros
        "(48726,'sf',0000978,0.3,1,'','',NULL,NULL);\n"
    )
    with gzip.open(path, "wt") as fh:
        fh.write(sql)
    return path


@pytest.fixture
def flat(tmp_path):
    path = tmp_path / "Domain2GO_supported_only_by_all.txt"
    path.write_text(
        "#Created at 3/04/2016 14:38:13\n"
        "Domain_type\tDomain_sunid\tGO_id\tGO_name\tGO_subontologies\t"
        "Information_content\tAnnotation_origin (Direct:1, Inherited:0)\n"
        "sf\t52540\tGO:0005524\tATP binding\tmolecular_function\t2.1\t1\n"
        "sf\t52540\tGO:0003674\tmolecular_function\tmolecular_function\t0.0\t0\n"
        "fa\t54870\tGO:0009055\telectron transfer\tmolecular_function\t3.0\t1\n"
    )
    return path


@pytest.fixture
def sp2go(tmp_path):
    path = tmp_path / "SP2GO.txt"
    path.write_text(
        "#Created at 3/04/2016 21:55:22\n"
        "Type\tSupra (comma-seperated sunids)\tGO_id\tGO_name\tGO_subontologies\t"
        "Information_content\tAnnotation_origin (Direct:1, Inherited:0)\n"
        "supra_sf\t46565,57938\tGO:0055035\tplastid\tcellular_component\t2.6\t1\n"
        "supra_sf\t46565,57938\tGO:0009536\tplastid\tcellular_component\t1.9\t0\n"
        # A single-component "supra" row: really a single domain, must be dropped
        # or the single-domain table would be double-counted.
        "supra_sf\t57802\tGO:0055035\tplastid\tcellular_component\t2.6\t1\n"
        "supra_sf\t46565,57938,49493\tGO:0016020\tmembrane\tcellular_component\t1.0\t1\n"
    )
    return path


class TestParseGoMapping:
    def test_keeps_superfamily_rows_only(self, dump):
        fdr, _hscore, sf_sunids, fa_sunids = compare.parse_go_mapping(dump)

        assert sf_sunids == {52540, 48726}
        assert fa_sunids == {54870}
        assert all(sunid in sf_sunids for sunid, _go in fdr)

    def test_reconstructs_zerofilled_go_ids(self, dump):
        fdr, _hscore, _sf, _fa = compare.parse_go_mapping(dump)

        assert (48726, "GO:0000978") in fdr
        assert (52540, "GO:0003674") in fdr

    def test_uses_all_score_not_single_score(self, dump):
        """``all_score`` is the comparator — we impose no single-domain filter."""
        fdr, _hscore, _sf, _fa = compare.parse_go_mapping(dump)

        assert fdr[(52540, "GO:0003674")] == pytest.approx(1e-8)  # not 0.5
        assert fdr[(52540, "GO:0005524")] == pytest.approx(0.002)  # not 1

    def test_hscore_only_where_not_null(self, dump):
        _fdr, hscore, _sf, _fa = compare.parse_go_mapping(dump)

        assert hscore == {(52540, "GO:0005524"): pytest.approx(95.45)}

    def test_commas_inside_quoted_text_do_not_shift_columns(self, dump):
        """``inherited_from_all`` holds comma-joined ids; a naive split breaks."""
        fdr, hscore, _sf, _fa = compare.parse_go_mapping(dump)

        # If the regex had mis-parsed the quoted "0042995,0009288" field, the
        # h-score would have been read off the wrong column (or missed).
        assert hscore[(52540, "GO:0005524")] == pytest.approx(95.45)
        assert fdr[(52540, "GO:0005524")] == pytest.approx(0.002)


class TestParseDomain2GoFlat:
    def test_splits_direct_from_inherited(self, flat):
        direct, inherited = compare.parse_domain2go_flat(flat)

        assert direct == {(52540, "GO:0005524")}
        assert inherited == {(52540, "GO:0003674")}

    def test_family_rows_excluded(self, flat):
        direct, inherited = compare.parse_domain2go_flat(flat)

        assert all(sunid != 54870 for sunid, _go in direct | inherited)

    def test_header_line_is_not_a_row(self, flat):
        direct, inherited = compare.parse_domain2go_flat(flat)

        assert len(direct) + len(inherited) == 2


class TestParseSp2go:
    def test_drops_single_component_rows(self, sp2go):
        supra = compare.parse_sp2go(sp2go, direct_only=False)

        assert ("57802",) not in supra
        assert set(supra) == {("46565", "57938"), ("46565", "57938", "49493")}

    def test_direct_only_filter(self, sp2go):
        direct = compare.parse_sp2go(sp2go, direct_only=True)
        both = compare.parse_sp2go(sp2go, direct_only=False)

        assert direct[("46565", "57938")] == {"GO:0055035"}
        assert both[("46565", "57938")] == {"GO:0055035", "GO:0009536"}


class TestGoNamespaces:
    def test_namespace_alt_id_and_obsolete(self, tmp_path):
        obo = tmp_path / "go-basic.obo"
        obo.write_text(
            "format-version: 1.2\n\n"
            "[Term]\nid: GO:0000001\nname: a\nnamespace: biological_process\n"
            "alt_id: GO:0099999\n\n"
            "[Term]\nid: GO:0000002\nname: b\nnamespace: molecular_function\n"
            "is_obsolete: true\n\n"
            "[Typedef]\nid: part_of\n"
        )
        namespaces, obsolete = compare.parse_go_namespaces(obo)

        assert namespaces["GO:0000001"] == "biological_process"
        # alt_ids must resolve, or a term merged since 2016 looks like a miss.
        assert namespaces["GO:0099999"] == "biological_process"
        assert obsolete == {"GO:0000002"}


class TestFisherPairs:
    def test_matches_a_hand_computed_table(self):
        domain_proteins = {"SSF1": {"P1", "P2", "P3"}}
        term_proteins = {"GO:1": {"P1", "P2", "P4"}}
        pvalues, overlaps = compare.fisher_p_for_pairs(
            [("SSF1", "GO:1")], domain_proteins, term_proteins, 10, n_jobs=1
        )

        assert overlaps[0] == 2
        # a=2 b=1 c=1 d=6; right-tailed Fisher on that table.
        from scipy.stats import fisher_exact

        expected = fisher_exact([[2, 1], [1, 6]], alternative="greater")[1]
        assert pvalues[0] == pytest.approx(expected, rel=1e-9)

    def test_absent_domain_or_term_gives_zero_overlap(self):
        pvalues, overlaps = compare.fisher_p_for_pairs(
            [("SSF404", "GO:404")], {}, {}, 10, n_jobs=1
        )

        assert overlaps[0] == 0
        assert pvalues[0] == pytest.approx(1.0)


class TestAgreement:
    def test_precision_recall_jaccard(self):
        shared, precision, recall, jaccard = compare.agreement({1, 2, 3}, {2, 3, 4, 5})

        assert shared == 2
        assert precision == "0.6667"
        assert recall == "0.5000"
        assert jaccard == "0.4000"

    def test_empty_sides_do_not_divide_by_zero(self):
        assert compare.agreement(set(), {1}) == (0, "", "0.0000", "0.0000")
