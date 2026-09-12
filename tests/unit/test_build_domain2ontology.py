"""Unit tests for scripts/build_domain2ontology.py (the release annotator).

The parsers are exercised on tiny synthetic fixtures that reproduce the
format quirks the real sources exhibit: wrapped DE/CC lines, alt_ids, EFO's
idspace-prefixed ids, ComplexPortal isoform tags, TCDB family inheritance,
Rhea directional ids, and the composite comma-joined term fallback.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPT_PATH = (
    Path(__file__).resolve().parents[2] / "scripts" / "build_domain2ontology.py"
)


@pytest.fixture(scope="module")
def build():
    spec = importlib.util.spec_from_file_location("build_domain2ontology", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules["build_domain2ontology"] = module
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# OBO parsing
# ---------------------------------------------------------------------------


OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: mitochondrion inheritance
alt_id: GO:0009999

[Term]
id: GO:0000002
name: obsolete mitochondrial genome maintenance
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""


def test_obo_names_primary_alt_and_obsolete(build, tmp_path):
    path = tmp_path / "mini.obo"
    path.write_text(OBO)
    names = build.parse_obo_names(path)
    assert names["GO:0000001"] == "mitochondrion inheritance"
    assert names["GO:0009999"] == "mitochondrion inheritance"  # alt_id folded in
    assert names["GO:0000002"].startswith("obsolete ")  # kept verbatim
    assert "part_of" not in names  # Typedef stanzas ignored


def test_obo_names_primary_wins_over_alt(build, tmp_path):
    path = tmp_path / "clash.obo"
    path.write_text(
        "[Term]\nid: X:1\nname: one\nalt_id: X:2\n\n[Term]\nid: X:2\nname: two\n"
    )
    names = build.parse_obo_names(path)
    assert names["X:2"] == "two"  # primary definition beats the alt_id row


def test_obo_names_efo_normalisation_and_underscore_variants(build, tmp_path):
    path = tmp_path / "efo.obo"
    path.write_text(
        "[Term]\nid: efo:EFO_0000001\nname: experimental factor\n\n"
        "[Term]\nid: OBA:VT0000217\nname: blood glucose amount\n"
    )
    from src.gwas_annotation_source import normalise_efo_id

    names = build.parse_obo_names(
        path, normalise=normalise_efo_id, underscore_variants=True
    )
    assert names["EFO:0000001"] == "experimental factor"
    assert names["OBA:VT0000217"] == "blood glucose amount"
    assert names["OBA_VT0000217"] == "blood glucose amount"  # GWAS URI spelling


# ---------------------------------------------------------------------------
# EC: enzyme.dat + enzclass.txt
# ---------------------------------------------------------------------------


ENZYME_DAT = """ID   1.1.1.1
DE   Alcohol dehydrogenase.
//
ID   1.1.1.2
DE   Very long enzyme name that Expasy wraps over
DE   two DE lines.
//
ID   1.1.1.3
DE   Transferred entry: 1.1.1.1.
//
ID   1.13.11.25
DE   3,4-dihydroxy-9,10-secoandrosta-1,3,5(10)-triene-9,17-dione 4,5-
DE   dioxygenase.
//
"""

ENZCLASS = """1. -. -.-  Oxidoreductases.
1. 1. -.-   Acting on the CH-OH group of donors.
1. 1. 1.-    With NAD(+) or NADP(+) as acceptor.
"""


def test_ec_names(build, tmp_path):
    enzyme = tmp_path / "enzyme.dat"
    enzyme.write_text(ENZYME_DAT)
    enzclass = tmp_path / "enzclass.txt"
    enzclass.write_text(ENZCLASS)
    names = build.parse_ec_names(enzyme, enzclass)
    assert names["1.1.1.1"] == "Alcohol dehydrogenase"
    assert (
        names["1.1.1.2"] == "Very long enzyme name that Expasy wraps over two DE lines"
    )
    assert "1.1.1.3" not in names  # Transferred placeholder skipped
    # Expasy wraps mid-token after a hyphen — the rejoin must not add a space.
    assert names["1.13.11.25"] == (
        "3,4-dihydroxy-9,10-secoandrosta-1,3,5(10)-triene-9,17-dione 4,5-dioxygenase"
    )
    assert names["1.-.-.-"] == "Oxidoreductases"
    assert names["1.1.-.-"] == "Acting on the CH-OH group of donors"
    assert names["1.1.1.-"] == "With NAD(+) or NADP(+) as acceptor"


# ---------------------------------------------------------------------------
# The Swiss-Prot pass
# ---------------------------------------------------------------------------


SPROT = """ID   TEST1_HUMAN             Reviewed;         100 AA.
AC   P00001;
CC   -!- FUNCTION: Something else entirely.
CC   -!- DISEASE: Immunodeficiency 104, severe combined, T-cell-
CC       negative (IMD104) [MIM:608971]: A disease. {ECO:0000269}.
CC   -!- DISEASE: Note=This block is a note, not a disease name.
CC   -!- CATALYTIC ACTIVITY:
CC       Reaction=ATP + H2O = ADP + phosphate + H(+); Xref=Rhea:RHEA:13065,
CC       ChEBI:CHEBI:15377; EC=3.6.4.6;
CC       PhysiologicalDirection=left-to-right; Xref=Rhea:RHEA:13066;
CC   ---------------------------------------------------------------------
DR   ComplexPortal; CPX-1621; TDRD3-TOP3B type IA topoisomerase complex. [Q9H7E2-3].
DR   DrugBank; DB00001; Lepirudin.
DR   TCDB; 8.A.98.1.6; the 14-3-3 protein (14-3-3) family.
DR   Reactome; R-HSA-140834; Extrinsic Pathway of Fibrin Clot Formation.
DR   MEROPS; S01.010; -.
//
ID   TEST2_HUMAN             Reviewed;         100 AA.
AC   P00002;
CC   -!- DISEASE: Alzheimer disease 1 (AD1) [MIM:104300]: A
CC       neurodegenerative disorder. {ECO:0000269}.
CC   -!- CATALYTIC ACTIVITY:
CC       Reaction=UDP + O-phospho-L-seryl-
CC       [protein] = water; Xref=Rhea:RHEA:20800, ChEBI:CHEBI:1;
DR   TCDB; 8.A.98.1.7; the 14-3-3 protein (14-3-3) family.
DR   ComplexPortal; CPX-365; Insulin hexamer.
DR   Reactome; R-DME-210693; STAT92E dimer transported to the cytosol. [Q9W0G1-2].
//
"""


@pytest.fixture(scope="module")
def sprot_tables(build, tmp_path_factory):
    path = tmp_path_factory.mktemp("sprot") / "uniprot_sprot.dat.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(SPROT)
    return build.build_swissprot_tables(path)


def test_disease_names(build, sprot_tables):
    # Hyphen-aware re-flow: "T-cell-" + "negative" joins without a space;
    # the trailing "(IMD104)" abbreviation and the [MIM:...] tail are cut.
    assert sprot_tables["disease"]["608971"] == (
        "Immunodeficiency 104, severe combined, T-cell-negative"
    )
    assert sprot_tables["disease"]["104300"] == "Alzheimer disease 1"
    assert len(sprot_tables["disease"]) == 2  # the Note= block minted nothing


def test_rhea_names(build, sprot_tables):
    assert sprot_tables["rhea"]["RHEA:13065"] == "ATP + H2O = ADP + phosphate + H(+)"
    assert sprot_tables["rhea"]["RHEA:13066"] == (
        "ATP + H2O = ADP + phosphate + H(+) (left-to-right)"
    )
    # A wrap after a hyphen before a '[protein]'-style token rejoins clean.
    assert sprot_tables["rhea"]["RHEA:20800"] == (
        "UDP + O-phospho-L-seryl-[protein] = water"
    )


def test_dr_names(build, sprot_tables):
    # ComplexPortal isoform tag stripped; chemistry brackets untouched.
    assert sprot_tables["complex"]["CPX-1621"] == (
        "TDRD3-TOP3B type IA topoisomerase complex"
    )
    assert sprot_tables["complex"]["CPX-365"] == "Insulin hexamer"
    assert sprot_tables["drugbank"]["DB00001"] == "Lepirudin"
    assert sprot_tables["reactome_dr"]["R-HSA-140834"] == (
        "Extrinsic Pathway of Fibrin Clot Formation"
    )
    # Reactome DR lines strip the isoform tag like ComplexPortal ones.
    assert sprot_tables["reactome_dr"]["R-DME-210693"] == (
        "STAT92E dimer transported to the cytosol"
    )


def test_tcdb_family_inheritance(build, sprot_tables):
    tcdb = sprot_tables["tcdb"]
    family = "the 14-3-3 protein (14-3-3) family"
    assert tcdb["8.A.98.1.6"] == family  # exact DR id
    assert tcdb["8.A.98"] == family  # family prefix inherits
    assert tcdb["8.A.98.1"] == family  # subfamily prefix inherits
    assert "8.A" not in tcdb  # class/subclass stay unnamed
    assert "8" not in tcdb


# ---------------------------------------------------------------------------
# Small vocabularies
# ---------------------------------------------------------------------------


def test_oncotree_names(build, tmp_path):
    path = tmp_path / "oncotree.json"
    path.write_text(
        json.dumps([{"code": "MT", "name": "Malignant Tumor", "parent": "TISSUE"}])
    )
    names = build.parse_oncotree_names(path)
    assert names["MT"] == "Malignant Tumor"
    assert names["TISSUE"] == "Tissue (OncoTree root)"  # virtual root covered


def test_cazy_names(build):
    assert build.cazy_name("GH27") == "Glycoside Hydrolase family GH27"
    assert build.cazy_name("GH") == "Glycoside Hydrolase family"  # interior node
    assert build.cazy_name("XY99") is None  # unknown prefix falls through


def test_merops_and_identity_and_nameless_providers(build, tmp_path):
    sources = build.NameSources(inputs={}, cache_dir=tmp_path)
    assert sources.provider("merops")("S") == "Serine peptidase"
    assert sources.provider("merops")("S01.010") is None  # family/leaf unnamed
    assert sources.provider("keyword")("Blood coagulation") == "Blood coagulation"
    assert sources.provider("unipathway")("UPA00109") is None
    with pytest.raises(KeyError):
        sources.provider("nonexistent")


def test_composite_term_fallback(build):
    table = {"WBbt:0004799": "F cell", "WBbt:0004942": "U cell"}
    assert build.named(table.get, "WBbt:0004799,WBbt:0004942") == "F cell / U cell"
    assert build.named(table.get, "WBbt:0004799,WBbt:missing") is None


def test_orphanet_names(build, tmp_path):
    owl = tmp_path / "ordo.owl"
    owl.write_text(
        '<Class rdf:about="http://www.orpha.net/ORDO/Orphanet_90349">\n'
        "  <rdfs:label>Acquired hemophilia &amp; friends</rdfs:label>\n"
        "</Class>\n"
        '<Class rdf:about="http://www.orpha.net/ORDO/Orphanet_C016">\n'
        "  <rdfs:label>metadata property, must be excluded</rdfs:label>\n"
        "</Class>\n"
    )
    xml = tmp_path / "en_product6.xml"
    xml.write_text(
        "<Disorder><OrphaCode>715694</OrphaCode>\n"
        "<ExpertLink>http://example.org</ExpertLink>\n"
        '<Name lang="en">Infant-type hemispheric glioma</Name></Disorder>\n'
        "<Disorder><OrphaCode>90349</OrphaCode>\n"
        '<Name lang="en">Should not override ORDO</Name></Disorder>\n'
    )
    names = build.parse_orphanet_names(owl, xml)
    assert names["90349"] == "Acquired hemophilia & friends"  # ORDO wins, unescaped
    assert names["715694"] == "Infant-type hemispheric glioma"  # XML fills the gap
    assert "C016" not in names and not any("C016" in key for key in names)


# ---------------------------------------------------------------------------
# End-to-end on a synthetic cell
# ---------------------------------------------------------------------------


HEADER = (
    "domain\tdemo_term\tp_value\tadj_p_value\todds_ratio\todds_ratio_ci_low\t"
    "odds_ratio_ci_high\thyper_score\tdomain_type\tconstituent_domains\t"
    "n_observations\ta\tb\tc\td\tic"
)


def test_annotate_cell_end_to_end(build, tmp_path, monkeypatch):
    cell_dir = tmp_path / "production" / "demo_baseline"
    cell_dir.mkdir(parents=True)
    stats = "\t".join(["1e-5", "1e-3", "2.0", "1.0", "4.0", "50.0"])
    tail = "\t".join(["3", "3", "1", "5", "100", "2.5"])
    (cell_dir / "domain_demo_associations_significant.tsv").write_text(
        f"{HEADER}\n"
        f"IPR000001\tX:1\t{stats}\tsingle\t-\t{tail}\n"
        f"IPR000001,IPR000002\tX:2\t{stats}\tsupra_pair\tIPR000001,IPR000002\t{tail}\n"
    )
    manifest = {
        "parameters": {
            "ontology": "keyword",
            "species": "human",
            "evidence_filter": "manual",
            "propagate_annotations": False,
        },
        "outputs": [
            {"role": "significant_associations", "sha256": "ab" * 32},
        ],
        "git": {"commit": "deadbeef"},
    }
    (cell_dir / "run_manifest_demo.json").write_text(json.dumps(manifest))

    out_dir = tmp_path / "final"
    out_dir.mkdir()
    sources = build.NameSources(inputs={}, cache_dir=out_dir / "names")
    domain_names = {"IPR000001": "Kringle", "IPR000002": "Serpin family"}
    result = build.annotate_cell(
        "demo_baseline",
        cell_dir / "domain_demo_associations_significant.tsv",
        manifest,
        sources,
        domain_names,
        out_dir,
    )

    with gzip.open(out_dir / "domain2-demo_baseline.tsv.gz", "rt") as handle:
        lines = [line.rstrip("\n").split("\t") for line in handle]
    assert lines[0][:4] == ["domain", "domain_name", "demo_term", "term_name"]
    assert lines[0][4:] == HEADER.split("\t")[2:]  # stat columns untouched
    assert lines[1][:4] == ["IPR000001", "Kringle", "X:1", "X:1"]  # keyword=identity
    assert lines[2][1] == "Kringle + Serpin family"  # supra join
    assert lines[1][4:] == [*stats.split("\t"), "single", "-", *tail.split("\t")]

    assert result.rows == 2
    assert result.distinct_domains == 2
    assert result.distinct_terms == 2
    assert result.named_terms == 2
    assert result.config == "baseline"
    assert result.sha256 == "ab" * 32


def test_annotate_cell_unnamed_terms_left_empty(build, tmp_path):
    cell_dir = tmp_path / "production" / "unipathway_baseline"
    cell_dir.mkdir(parents=True)
    stats = "\t".join(["1e-5", "1e-3", "2.0", "1.0", "4.0", "50.0"])
    tail = "\t".join(["3", "3", "1", "5", "100", "2.5"])
    (cell_dir / "domain_unipathway_associations_significant.tsv").write_text(
        f"{HEADER}\nIPR000001\tUPA00109\t{stats}\tsingle\t-\t{tail}\n"
    )
    manifest = {"parameters": {"ontology": "unipathway"}, "outputs": []}
    out_dir = tmp_path / "final"
    out_dir.mkdir()
    sources = build.NameSources(inputs={}, cache_dir=out_dir / "names")
    result = build.annotate_cell(
        "unipathway_baseline",
        cell_dir / "domain_unipathway_associations_significant.tsv",
        manifest,
        sources,
        {},
        out_dir,
    )
    with gzip.open(out_dir / "domain2-unipathway_baseline.tsv.gz", "rt") as handle:
        lines = [line.rstrip("\n").split("\t") for line in handle]
    assert lines[1][3] == ""  # unnamed term → empty column, row kept
    assert lines[1][1] == "IPR000001"  # unknown domain falls back to its id
    assert result.named_terms == 0


def test_name_cache_roundtrip(build, tmp_path):
    # Whitespace normalisation happens at build time so a cold (building) run
    # and a warm (cache-reading) run annotate with byte-identical names.
    sources = build.NameSources(inputs={}, cache_dir=tmp_path / "names")
    table = sources._cached("demo", lambda: {"A:1": "alpha  beta\tgamma"})
    assert table["A:1"] == "alpha beta gamma"
    # A fresh instance must read the cache, not rebuild (builder would raise).
    reread = build.NameSources(inputs={}, cache_dir=tmp_path / "names")
    cached = reread._cached("demo", lambda: pytest.fail("cache not used"))
    assert cached["A:1"] == table["A:1"]


SUBCELL = """ID   Cell membrane.
AC   SL-0039
//
IT   Multi-pass membrane protein.
AC   SL-9909
//
IO   Extracellular side.
AC   SL-9906
//
"""


def test_subcellular_names(build, tmp_path):
    path = tmp_path / "subcell.txt"
    path.write_text(SUBCELL)
    names = build.parse_subcellular_names(path)
    assert names["SL-0039"] == "Cell membrane"
    assert names["SL-9909"] == "Multi-pass membrane protein"  # IT entries share
    assert names["SL-9906"] == "Extracellular side"  # the SL space, as do IO
