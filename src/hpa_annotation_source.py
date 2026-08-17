"""HPA single-cell cell-type expression: elevated genes, re-keyed to UniProt.

The Human Protein Atlas ``rna_single_cell_type.tsv`` gives one expression value
(nCPM) per gene × cell type over 154 single-cell cell types. This is an
**expression** layer: an annotation means "elevated in this cell type", not
"phenotype of" — the loosest semantics after the GWAS layer, stated here so a
domain → cell-type association is read as enrichment among markers of that
cell type, nothing more.

Elevation policy
----------------
"Expressed at all" is not an annotation: at the HPA detection cutoff
(nCPM ≥ 1) 1,940,503 of the 3,087,080 gene × cell-type pairs qualify (63% —
nearly every gene "annotated" to most cell types, at acquisition, HPA
2025-12-12). The layer therefore keeps only *elevated* expression, using HPA's
own "cell type enhanced" criterion: ``nCPM ≥ 1`` **and** ``nCPM ≥ 4× the mean
of the gene's other cell types``. (HPA's "enriched" class — ≥ 4× every other
cell type — is a subset, so both of HPA's elevated classes pass; "low
specificity" does not.) That keeps 97,023 pairs over 17,407 genes — ~5.6 cell
types per gene, a usable annotation density. Both the raw and kept counts are
returned (:class:`HPAFilterCounts`) and logged.

Cell Ontology assessment (why the terms are HPA names, not CL ids)
------------------------------------------------------------------
The brief for this layer was CL. The HPA file carries no CL ids, and no
name→CL mapping ships in the download, so the only honest route was matching
HPA cell-type names against ``cl.obo`` names and EXACT synonyms — exact and
case-insensitive matching only (plus trivial trailing-``s`` singularisation),
no fuzzy matching. Result against CL 2026-06-08: **73 of 154 cell types
(47%)**, below the 60% floor set for delivering a CL-keyed layer. The misses
are systematic, not typographic — HPA qualifies names by tissue ("salivary
duct cells", "endometrial ciliated cells") where CL's classes are unqualified
or split differently, and mapping those correctly is curation, not string
matching. The layer is therefore delivered **flat**, keyed by HPA's own
cell-type vocabulary (no hierarchy, so no True Path or relative inference);
re-keying onto CL stays open until a curated mapping exists.

The gene axis is Ensembl gene ids, re-keyed to UniProt accessions at parse
time via the flat file's ``DR Ensembl`` gene field
(:func:`src.gene_mapping.parse_ensembl_gene_accession_map`), with the counted
policy for unmapped and one-to-many ids.
"""

from __future__ import annotations

import csv
import io
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec
from src.gene_mapping import (
    parse_ensembl_gene_accession_map,
    remap_gene_annotations,
)
from src.remap import RemapCoverage

#: HPA single-cell cell types, by name (flat — see the module docstring for
#: why these are not CL ids).
CELLTYPE_SPEC = OntologySpec(
    ontology_id="HPA-CellType", name="Human Protein Atlas single-cell cell type"
)

#: HPA's detection cutoff: below this a gene is "not expressed" in the type.
MIN_EXPRESSION = 1.0

#: HPA's "cell type enhanced" fold: elevated means ≥ this × the mean of the
#: gene's other cell types (their "enriched" class is a subset of this).
ELEVATION_FOLD = 4.0

#: HPA has renamed the expression column between releases.
_EXPRESSION_COLUMNS = ("nCPM", "nTPM")


@dataclass(frozen=True)
class HPAFilterCounts:
    """What the elevation policy did to the expression matrix.

    Attributes:
        n_pairs: gene × cell-type rows read.
        n_expressed: pairs at or above the detection cutoff (``nCPM ≥ 1``) —
            what "expressed in" would have annotated, kept as the honest
            comparison point for the elevation filter.
        n_elevated: pairs passing the elevation criterion (the annotations).
        n_genes: distinct genes read.
        n_genes_elevated: genes with at least one elevated cell type.
    """

    n_pairs: int = 0
    n_expressed: int = 0
    n_elevated: int = 0
    n_genes: int = 0
    n_genes_elevated: int = 0


def _open_expression_table(path: Path):
    """Open the expression TSV, reaching inside the distribution zip if needed."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"HPA single-cell expression file not found: {path}")
    if path.suffix == ".zip":
        archive = zipfile.ZipFile(path)
        members = archive.namelist()
        if len(members) != 1:
            raise ValueError(
                f"{path} should contain exactly the expression TSV; found {members}"
            )
        return io.TextIOWrapper(archive.open(members[0]), encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def parse_hpa_single_cell(
    path: Path,
) -> Tuple[Dict[str, Set[str]], HPAFilterCounts]:
    """Parse HPA's single-cell table into ``{Ensembl gene: {cell type}}``.

    Keeps a (gene, cell type) pair only under the elevation policy in the
    module docstring. The expression column is found by name (``nCPM``, or the
    older ``nTPM``); unparsable values are treated as 0 (HPA writes ``0.0``
    for absence, so nothing legitimate is lost).
    """
    logger.info(f"Parsing HPA single-cell expression from {path}")
    per_gene: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
    with _open_expression_table(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError(f"{path} is empty")
        try:
            gene_col = header.index("Gene")
            type_col = header.index("Cell type")
        except ValueError as exc:
            raise ValueError(f"{path} header {header} lacks {exc}") from exc
        expression_col = next(
            (header.index(name) for name in _EXPRESSION_COLUMNS if name in header),
            None,
        )
        if expression_col is None:
            raise ValueError(
                f"{path} header {header} has none of the expression columns "
                f"{_EXPRESSION_COLUMNS}"
            )
        for row in reader:
            if len(row) <= max(gene_col, type_col, expression_col):
                continue
            gene = row[gene_col].strip()
            cell_type = row[type_col].strip()
            if not gene or not cell_type:
                continue
            try:
                value = float(row[expression_col])
            except ValueError:
                value = 0.0
            per_gene[gene].append((cell_type, value))

    gene_terms: Dict[str, Set[str]] = {}
    n_pairs = n_expressed = n_elevated = 0
    for gene, values in per_gene.items():
        n_pairs += len(values)
        total = sum(value for _, value in values)
        elevated: Set[str] = set()
        for cell_type, value in values:
            if value >= MIN_EXPRESSION:
                n_expressed += 1
                others_mean = (
                    (total - value) / (len(values) - 1) if len(values) > 1 else 0.0
                )
                if value >= ELEVATION_FOLD * others_mean:
                    elevated.add(cell_type)
        if elevated:
            n_elevated += len(elevated)
            gene_terms[gene] = elevated

    counts = HPAFilterCounts(
        n_pairs=n_pairs,
        n_expressed=n_expressed,
        n_elevated=n_elevated,
        n_genes=len(per_gene),
        n_genes_elevated=len(gene_terms),
    )
    logger.info(
        f"  Pairs: {counts.n_pairs:,}; expressed (≥ {MIN_EXPRESSION:g}): "
        f"{counts.n_expressed:,}; elevated (≥ {ELEVATION_FOLD:g}× mean of other "
        f"cell types): {counts.n_elevated:,} over "
        f"{counts.n_genes_elevated:,} / {counts.n_genes:,} genes"
    )
    return gene_terms, counts


class HPACellTypeAnnotationSource(AnnotationSource):
    """Domain annotations keyed by HPA single-cell cell-type name.

    Reads the expression table under the elevation policy and re-keys its
    Ensembl gene ids to UniProt accessions before the statistics see them.
    """

    spec = CELLTYPE_SPEC

    def __init__(self, expression_path: Path, dat_path: Path) -> None:
        self.expression_path = Path(expression_path)
        self.dat_path = Path(dat_path)
        #: Populated by :meth:`parse`: the elevation-policy audit.
        self.filter_counts: Optional[HPAFilterCounts] = None
        #: Populated by :meth:`parse`; its *values* range over Ensembl gene
        #: ids, its *keys* over cell-type names (axis-swapped, see
        #: remap_gene_annotations).
        self.coverage: Optional[RemapCoverage] = None

    def parse(self) -> Dict[str, Set[str]]:
        gene_terms, self.filter_counts = parse_hpa_single_cell(self.expression_path)
        gene_map = parse_ensembl_gene_accession_map(self.dat_path)
        remapped, self.coverage = remap_gene_annotations(
            gene_terms, gene_map, label="Ensembl gene→UniProt (HPA cell type)"
        )
        return remapped
