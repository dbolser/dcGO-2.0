"""Human Phenotype Ontology: gene-level phenotype annotations, re-keyed to UniProt.

HPO distributes ``genes_to_phenotype.txt`` — a TSV of NCBI GeneID → HP term,
derived from the disease-level ``phenotype.hpoa`` via each disease's gene
associations. Its keys are genes, not proteins, so the layer re-keys
GeneID → UniProt accession *at parse time* using the ``DR   GeneID`` lines of
the Swiss-Prot flat file (:mod:`src.gene_mapping`) — the same
translate-before-the-statistics pattern as the DOID layer, with the same
counted policy for unmapped and one-to-many ids.

The hierarchy is ``hp.obo``, read by the shared light OBO reader
(:func:`src.hierarchy.parse_obo_child_parents`); no subtree restriction is
applied — dcGO layers are deliberately unbiased, and inheritance-mode or
clinical-modifier terms simply behave like any other sparse term.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Set

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec
from src.disease_ontology import RemapCoverage
from src.gene_mapping import parse_gene_accession_index, remap_gene_annotations

#: Human Phenotype Ontology terms, e.g. ``HP:0001250``.
HPO_SPEC = OntologySpec(
    ontology_id="HP", name="Human Phenotype Ontology", term_prefix="HP:"
)


def parse_genes_to_phenotype(path: Path) -> Dict[str, Set[str]]:
    """Parse HPO's ``genes_to_phenotype.txt`` into ``{ncbi gene id: {HP term}}``.

    Columns are read by header name (``ncbi_gene_id``, ``hpo_id``), so the
    frequency/disease columns HPO adds or reorders between releases are
    ignored rather than positional-index hazards.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"HPO genes_to_phenotype file not found: {path}")

    logger.info(f"Parsing HPO gene→phenotype annotations from {path}")
    gene_terms: Dict[str, Set[str]] = defaultdict(set)
    n_rows = 0
    with open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"ncbi_gene_id", "hpo_id"} - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"{path} lacks expected column(s) {sorted(missing)}; "
                f"header was {reader.fieldnames}"
            )
        for row in reader:
            gene_id = (row["ncbi_gene_id"] or "").strip()
            hpo_id = (row["hpo_id"] or "").strip()
            if gene_id and hpo_id:
                n_rows += 1
                gene_terms[gene_id].add(hpo_id)
    logger.info(
        f"  Rows: {n_rows:,}; genes: {len(gene_terms):,}; terms: "
        f"{len({term for terms in gene_terms.values() for term in terms}):,}"
    )
    return dict(gene_terms)


class HPOAnnotationSource(AnnotationSource):
    """Domain annotations keyed by HPO term.

    Reads ``genes_to_phenotype.txt`` and re-keys its NCBI GeneIDs to UniProt
    accessions before the annotations reach the statistics.
    """

    spec = HPO_SPEC

    def __init__(self, genes_to_phenotype_path: Path, dat_path: Path) -> None:
        self.genes_to_phenotype_path = Path(genes_to_phenotype_path)
        self.dat_path = Path(dat_path)
        #: Populated by :meth:`parse`; axis-swapped as documented in
        #: :func:`src.gene_mapping.remap_gene_annotations`.
        self.coverage: Optional[RemapCoverage] = None

    def parse(self) -> Dict[str, Set[str]]:
        gene_terms = parse_genes_to_phenotype(self.genes_to_phenotype_path)
        index = parse_gene_accession_index(self.dat_path)
        remapped, self.coverage = remap_gene_annotations(
            gene_terms, index.geneid, label="GeneID→UniProt (HPO)"
        )
        return remapped
