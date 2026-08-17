"""Unit tests for the SynGO layer (HGNC genes re-keyed to UniProt).

The fixture builds a miniature bulk-release zip — real xlsx sheets, since that
is the only format SynGO ships — with one gene per mapping case: HGNC-id
mapped, symbol-fallback mapped, and unmapped.
"""

import zipfile

import openpyxl
import pytest

from src.gene_mapping import parse_gene_accession_index
from src.hierarchy import closure_ancestors
from src.syngo_annotation_source import (
    SynGOAnnotationSource,
    parse_syngo_annotations,
    parse_syngo_hierarchy,
    resolve_hgnc_accessions,
)

ANNOTATION_ROWS = [
    # (hgnc_id, hgnc_symbol, go_id). HGNC:1 maps by id; HGNC:2 only by symbol;
    # HGNC:9 maps to nothing. The uniprot_id column carries the (possibly
    # non-human) evidence protein and must be ignored.
    ("HGNC:1", "GENE1", "Q00000", "GO:0045202"),
    ("HGNC:1", "GENE1", "Q00000", "SYNGO:presynapse_x"),
    ("HGNC:2", "GENE2", "Q11111", "GO:0045202"),
    ("HGNC:9", "GENE9", "Q22222", "GO:0098793"),
    # A multi-gene row: a paralog pair the evidence could not separate.
    ("HGNC:3;HGNC:4", "GENE3;GENE4", "Q33333", "GO:0098793"),
]

ONTOLOGY_ROWS = [
    # (id, parent_id): a three-level chain under a parentless root.
    ("GO:0045202", None),
    ("GO:0098793", "GO:0045202"),
    ("SYNGO:presynapse_x", "GO:0098793"),
]

DAT = """AC   P20001;
DR   HGNC; HGNC:1; GENE1.
//
AC   P20002;
DR   HGNC; HGNC:22; GENE2.
//
"""


def _sheet(rows, header):
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.append(header)
    for row in rows:
        sheet.append(row)
    return workbook


@pytest.fixture
def zip_path(tmp_path):
    annotations = tmp_path / "annotations.xlsx"
    _sheet(ANNOTATION_ROWS, ["hgnc_id", "hgnc_symbol", "uniprot_id", "go_id"]).save(
        annotations
    )
    ontologies = tmp_path / "ontologies.xlsx"
    _sheet(ONTOLOGY_ROWS, ["id", "parent_id"]).save(ontologies)
    path = tmp_path / "syngo_complete_data.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.write(annotations, "annotations.xlsx")
        archive.write(ontologies, "ontologies.xlsx")
    return path


@pytest.fixture
def dat(tmp_path):
    path = tmp_path / "uniprot_sprot.dat"
    path.write_text(DAT)
    return path


class TestParseSynGOAnnotations:
    def test_terms_are_grouped_by_hgnc_id(self, zip_path):
        gene_terms, _ = parse_syngo_annotations(zip_path)
        assert gene_terms == {
            "HGNC:1": {"GO:0045202", "SYNGO:presynapse_x"},
            "HGNC:2": {"GO:0045202"},
            "HGNC:9": {"GO:0098793"},
            "HGNC:3": {"GO:0098793"},
            "HGNC:4": {"GO:0098793"},
        }

    def test_multi_gene_rows_credit_every_listed_gene(self, zip_path):
        gene_terms, symbols = parse_syngo_annotations(zip_path)
        assert gene_terms["HGNC:3"] == gene_terms["HGNC:4"] == {"GO:0098793"}
        assert (symbols["HGNC:3"], symbols["HGNC:4"]) == ("GENE3", "GENE4")

    def test_symbols_are_kept_for_the_fallback(self, zip_path):
        _, symbols = parse_syngo_annotations(zip_path)
        assert symbols == {
            "HGNC:1": "GENE1",
            "HGNC:2": "GENE2",
            "HGNC:9": "GENE9",
            "HGNC:3": "GENE3",
            "HGNC:4": "GENE4",
        }

    def test_missing_zip_fails_loudly(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_syngo_annotations(tmp_path / "gone.zip")


class TestParseSynGOHierarchy:
    def test_child_to_parent_edges(self, zip_path):
        assert parse_syngo_hierarchy(zip_path) == {
            "GO:0098793": {"GO:0045202"},
            "SYNGO:presynapse_x": {"GO:0098793"},
        }

    def test_closure_walks_to_the_root(self, zip_path):
        ancestors = closure_ancestors(parse_syngo_hierarchy(zip_path))
        assert ancestors("SYNGO:presynapse_x") == {"GO:0098793", "GO:0045202"}


class TestResolveHGNCAccessions:
    def test_symbol_fallback_is_used_and_counted(self, zip_path, dat):
        gene_terms, symbols = parse_syngo_annotations(zip_path)
        index = parse_gene_accession_index(dat)
        gene_map, n_fallback = resolve_hgnc_accessions(gene_terms, symbols, index)
        assert gene_map.targets("HGNC:1") == {"P20001"}  # by HGNC id
        assert gene_map.targets("HGNC:2") == {"P20002"}  # by symbol only
        assert gene_map.targets("HGNC:9") == set()
        assert n_fallback == 1


class TestSynGOAnnotationSource:
    """End-to-end: HGNC-keyed xlsx in, accession-keyed SynGO terms out."""

    def test_parse_produces_accession_keyed_terms(self, zip_path, dat):
        source = SynGOAnnotationSource(zip_path, dat)
        assert source.parse() == {
            "P20001": {"GO:0045202", "SYNGO:presynapse_x"},
            "P20002": {"GO:0045202"},
        }

    def test_unmapped_gene_and_fallback_are_counted(self, zip_path, dat):
        source = SynGOAnnotationSource(zip_path, dat)
        source.parse()
        assert source.coverage.unmapped_terms == ["HGNC:3", "HGNC:4", "HGNC:9"]
        assert source.n_symbol_fallback == 1

    def test_spec_has_no_single_term_prefix(self, zip_path, dat):
        # SynGO terms mix GO: and SYNGO: ids, so the spec cannot pin one.
        spec = SynGOAnnotationSource(zip_path, dat).spec
        assert spec.ontology_id == "SynGO"
        assert spec.term_prefix is None
