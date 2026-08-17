"""FlyBase allele phenotypes as FBcv/FBbt layers, re-keyed to UniProt.

The dcGO trick applied to *Drosophila* (``--species fly --ontology fbcv`` /
``fbbt``): domain → phenotype-term associations are learned on fly proteins,
and the species-agnostic domains carry them anywhere.

Annotation input is FlyBase's ``genotype_phenotype_data_*.tsv.gz``: one row per
genotype × phenotype term, columns (named in the ``#genotype_symbols`` header)
``genotype_symbols, genotype_FBids, phenotype_name, phenotype_id,
qualifier_names, qualifier_ids, reference``. The phenotype column mixes two
vocabularies — FBcv phenotype classes ("abnormal learning") and FBbt anatomy
terms (the structure manifesting the phenotype) — plus a minority of GO/SO
ids. One parser serves both registry keys; each keeps only its own prefix and
counts the rest as skipped. The FBbt layer's semantics are therefore
"abnormality manifested in this anatomical structure", the fly analogue of the
ZFA layer (:mod:`src.zfin_annotation_source`).

**Single-allele policy.** ``genotype_FBids`` lists every zygosity component,
``/``-separated within a locus and space-separated across loci
(``FBal0059629/FBal0408868 FBal0134436``). A phenotype observed on a
multi-allele genotype cannot be attributed to one allele, so only genotypes
whose id set is exactly one distinct FBal id are kept (a homozygote like
``FBal0119724`` counts once). At acquisition this keeps 199,987 of 399,972
data rows (50.0%); multi-allele rows and rows keyed by a non-allele id
(aberrations, insertions) are dropped and counted.

Kept rows are re-keyed FBal → FBgn via ``fbal_to_fbgn_*.tsv.gz`` (allele →
gene is 1:1 in FlyBase; misses are counted), pooling alleles of the same gene,
then FBgn → UniProt via the accession column of ``fbgn_NAseq_Uniprot_*.tsv.gz``
with the shared counted policy for unmapped and one-to-many ids
(:func:`src.gene_mapping.remap_gene_annotations`). Both coverages are exposed.

Hierarchies are ``fbcv.obo`` / ``fbbt.obo``, read by the shared light OBO
reader (:func:`src.hierarchy.parse_obo_child_parents`).
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Set

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec
from src.gene_mapping import GeneAccessionMap, remap_gene_annotations
from src.remap import RemapCoverage

#: FlyBase controlled vocabulary phenotype classes, e.g. ``FBcv:0000397``.
FBCV_SPEC = OntologySpec(
    ontology_id="FBcv", name="FlyBase phenotype class", term_prefix="FBcv:"
)
#: Drosophila anatomy terms manifesting a phenotype, e.g. ``FBbt:00004230``.
FBBT_SPEC = OntologySpec(
    ontology_id="FBbt", name="Drosophila Anatomy Ontology", term_prefix="FBbt:"
)

_FBIDS_COL = 1
_TERM_COL = 3
_ID_SPLIT = re.compile(r"[ /]+")


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def parse_genotype_phenotype(path: Path, term_prefix: str) -> Dict[str, Set[str]]:
    """Parse FlyBase genotype phenotypes into ``{FBal id: {term}}``.

    Keeps only single-allele genotypes and terms carrying ``term_prefix`` —
    see the module docstring for the policy and the counted drops.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FlyBase genotype_phenotype file not found: {path}")

    logger.info(f"Parsing FlyBase phenotype annotations from {path}")
    allele_terms: Dict[str, Set[str]] = defaultdict(set)
    n_rows = 0
    n_multi_allele = 0
    n_non_allele = 0
    n_other_prefix = 0
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= _TERM_COL:
                continue
            ids = {token for token in _ID_SPLIT.split(fields[_FBIDS_COL]) if token}
            if len(ids) != 1:
                n_multi_allele += 1
                continue
            (allele,) = ids
            if not allele.startswith("FBal"):
                n_non_allele += 1
                continue
            term = fields[_TERM_COL].strip()
            if not term.startswith(term_prefix):
                n_other_prefix += 1
                continue
            n_rows += 1
            allele_terms[allele].add(term)
    logger.info(
        f"  {term_prefix} rows kept: {n_rows:,} (single-allele); dropped: "
        f"{n_multi_allele:,} multi-allele, {n_non_allele:,} non-allele ids, "
        f"{n_other_prefix:,} other-vocabulary terms"
    )
    return dict(allele_terms)


def parse_fbal_to_fbgn(path: Path) -> Dict[str, str]:
    """Parse ``fbal_to_fbgn`` into ``{FBal id: FBgn id}`` (1:1 in FlyBase)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FlyBase fbal_to_fbgn file not found: {path}")

    mapping: Dict[str, str] = {}
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0].startswith("FBal"):
                mapping[fields[0]] = fields[2]
    logger.info(f"  FBal → FBgn mappings: {len(mapping):,}")
    return mapping


def parse_fbgn_uniprot(path: Path) -> GeneAccessionMap:
    """Parse ``fbgn_NAseq_Uniprot`` into an FBgn → accessions map.

    Column 3 is the primary FBgn, column 6 the UniProt accession (most rows
    are nucleotide cross-references with an empty accession field and are
    skipped). One gene collects every accession it is paired with — kept
    one-to-many for the counted expansion policy downstream.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"FlyBase fbgn_NAseq_Uniprot file not found: {path}")

    logger.info(f"Building FBgn → accession map from {path}")
    mapping: Dict[str, Set[str]] = defaultdict(set)
    n_rows = 0
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            gene, accession = fields[2].strip(), fields[5].strip()
            if gene.startswith("FBgn") and accession:
                n_rows += 1
                mapping[gene].add(accession)
    gene_map = GeneAccessionMap("FBgn", dict(mapping), n_rows)
    logger.info(
        f"  Genes with accessions: {len(gene_map):,} "
        f"({gene_map.n_one_to_many:,} one-to-many)"
    )
    return gene_map


class FlyBasePhenotypeAnnotationSource(AnnotationSource):
    """Domain annotations keyed by FBcv or FBbt phenotype term.

    Reads the genotype-phenotype table (single-allele genotypes, one term
    prefix), pools alleles to genes via ``fbal_to_fbgn``, and re-keys FBgn ids
    to UniProt accessions before the statistics see them.
    """

    def __init__(
        self,
        genotype_phenotype_path: Path,
        fbal_to_fbgn_path: Path,
        fbgn_uniprot_path: Path,
        spec: OntologySpec = FBCV_SPEC,
    ) -> None:
        self.genotype_phenotype_path = Path(genotype_phenotype_path)
        self.fbal_to_fbgn_path = Path(fbal_to_fbgn_path)
        self.fbgn_uniprot_path = Path(fbgn_uniprot_path)
        self.spec = spec
        #: FBal ids whose gene was unknown to fbal_to_fbgn (dropped, counted).
        self.n_unmapped_alleles: int = 0
        #: Populated by :meth:`parse`; its *values* range over FBgn ids, its
        #: *keys* over phenotype terms (axis-swapped, see remap_gene_annotations).
        self.coverage: Optional[RemapCoverage] = None

    def parse(self) -> Dict[str, Set[str]]:
        allele_terms = parse_genotype_phenotype(
            self.genotype_phenotype_path, self.spec.term_prefix
        )
        fbal_to_fbgn = parse_fbal_to_fbgn(self.fbal_to_fbgn_path)

        gene_terms: Dict[str, Set[str]] = defaultdict(set)
        self.n_unmapped_alleles = 0
        for allele, terms in allele_terms.items():
            gene = fbal_to_fbgn.get(allele)
            if gene is None:
                self.n_unmapped_alleles += 1
                continue
            gene_terms[gene] |= terms
        logger.info(
            f"  Alleles pooled onto {len(gene_terms):,} genes; "
            f"{self.n_unmapped_alleles:,} alleles had no FBgn and were dropped"
        )

        gene_map = parse_fbgn_uniprot(self.fbgn_uniprot_path)
        remapped, self.coverage = remap_gene_annotations(
            dict(gene_terms),
            gene_map,
            label=f"FBgn→UniProt ({self.spec.ontology_id})",
            target_label="UniProt accession",
        )
        return remapped
