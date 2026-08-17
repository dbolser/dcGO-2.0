"""EFO trait annotations from the GWAS Catalog, re-keyed to UniProt.

The NHGRI-EBI GWAS Catalog's ontology-annotated association file maps SNP
associations to Experimental Factor Ontology traits (``MAPPED_TRAIT_URI``) and
to the genes the SNPs fall in (``MAPPED_GENE``). Collapsing the SNP axis gives
a gene → trait map: "variation in this gene is associated with this trait" —
*genetic association* semantics, not curated function, stated here because it
is the loosest evidence type of any layer in the registry.

Row policy (every drop counted and logged, and returned to callers as
:class:`GWASFilterCounts`):

* **Genome-wide significance** — only rows with ``P-VALUE < 5e-8`` are kept,
  the field-standard threshold. The catalog includes sub-threshold secondary
  signals (182,812 of 1,188,619 rows at the 2026-08-02 release, unparsable
  p-values included); counting them as annotations would let suggestive hits
  carry the same weight as replicated ones.
* **Intergenic associations** — dropped. A SNP between two genes cannot be
  attributed to either without a linkage model this pipeline does not have.
  Two forms exist: rows flagged ``INTERGENIC == 1`` (412,503 at acquisition)
  and rows whose ``MAPPED_GENE`` is a flanking-gene pair written ``"A - B"``
  (2,325 more; the spaced dash is the flanking notation — gene symbols
  containing hyphens, like ``HLA-B``, are never written with spaces).
* **Gene column** — ``MAPPED_GENE`` (Ensembl-pipeline mapping of the SNP
  location), not ``REPORTED GENE(S)`` (free-text author claims, frequently
  missing, ``"NR"`` or marketing names). Multiple mapped genes (``", "`` or
  ``"; "`` separated) each receive the trait.
* **Traits** — every ``MAPPED_TRAIT_URI`` (comma-separated when a study maps
  to several) is normalised from URI to CURIE (``.../EFO_0004574`` →
  ``EFO:0004574``). Trait CURIEs outside EFO's own namespace (OBA, MONDO, HP,
  GO …) are kept: ``efo.obo`` imports those terms, so they resolve in the
  hierarchy like any native EFO id. Rows with no trait URI are dropped and
  counted.

The gene symbol → UniProt re-key uses the Swiss-Prot flat file's ``DR HGNC``
symbol index (:mod:`src.gene_mapping`), with the counted policy for unmapped
and one-to-many symbols. EFO is large and messy — the unmapped-symbol count is
the honest measure of this layer, not the row count.

The hierarchy is ``efo.obo``. Its own term ids carry an OBO idspace prefix
(``id: efo:EFO_0000001``) while imported terms are plain CURIEs; the loader
normalises both sides of every edge to the CURIE form the annotations use.
"""

from __future__ import annotations

import csv
import io
import re
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec
from src.gene_mapping import parse_gene_accession_index, remap_gene_annotations
from src.hierarchy import parse_obo_child_parents
from src.remap import RemapCoverage

#: GWAS Catalog traits, EFO-keyed (with imported OBA/MONDO/HP/GO ids kept).
EFO_SPEC = OntologySpec(
    ontology_id="EFO", name="Experimental Factor Ontology", term_prefix="EFO:"
)

#: The genome-wide significance threshold (Risch & Merikangas convention).
GENOME_WIDE_SIGNIFICANCE = 5e-8

_REQUIRED_COLUMNS = ("MAPPED_GENE", "MAPPED_TRAIT_URI", "P-VALUE", "INTERGENIC")

#: An OBO idspace prefix (lowercase) in front of an underscore CURIE,
#: e.g. the ``efo:`` in ``efo:EFO_0000001``.
_IDSPACE_RE = re.compile(r"^[a-z][a-z0-9]*:(?=[A-Za-z]+_\d)")
_UNDERSCORE_CURIE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]*_\d+$")


def normalise_efo_id(term: str) -> str:
    """Normalise an EFO/OBO identifier to CURIE form.

    ``efo:EFO_0000001`` and ``EFO_0000001`` become ``EFO:0000001``;
    already-CURIE ids (``MONDO:0000001``, ``OBA:0000015``) and full URLs pass
    through unchanged.
    """
    term = _IDSPACE_RE.sub("", term)
    if _UNDERSCORE_CURIE_RE.fullmatch(term):
        term = term.replace("_", ":", 1)
    return term


def trait_uri_to_curie(uri: str) -> str:
    """Turn a ``MAPPED_TRAIT_URI`` into the CURIE the hierarchy uses.

    ``http://www.ebi.ac.uk/efo/EFO_0004574`` → ``EFO:0004574``;
    ``http://purl.obolibrary.org/obo/MONDO_0005148`` → ``MONDO:0005148``.
    """
    return normalise_efo_id(uri.rstrip("/").rsplit("/", 1)[-1])


def parse_efo_child_parents(path: Path) -> Dict[str, set]:
    """``efo.obo`` as a child→parents map, ids normalised to CURIE form.

    Uses the shared light OBO reader (``is_a`` + ``part_of``), then rewrites
    both ends of every edge with :func:`normalise_efo_id` so EFO's
    idspace-prefixed native ids join the same graph as its imported CURIEs.
    """
    raw = parse_obo_child_parents(path)
    return {
        normalise_efo_id(child): {normalise_efo_id(parent) for parent in parents}
        for child, parents in raw.items()
    }


@dataclass(frozen=True)
class GWASFilterCounts:
    """What the row policy did to the association file.

    Attributes:
        n_rows: data rows read.
        n_below_significance: rows dropped by the ``P-VALUE < 5e-8`` filter
            (unparsable p-values included).
        n_intergenic: rows flagged ``INTERGENIC == 1``.
        n_flanking: further rows whose ``MAPPED_GENE`` was a flanking
            ``"A - B"`` pair (or empty) — intergenic in effect.
        n_no_trait: rows with no ``MAPPED_TRAIT_URI``.
        n_kept: rows contributing at least one (gene, trait) pair.
    """

    n_rows: int = 0
    n_below_significance: int = 0
    n_intergenic: int = 0
    n_flanking: int = 0
    n_no_trait: int = 0
    n_kept: int = 0


def _open_association_table(path: Path):
    """Open the association TSV, reaching inside the distribution zip if needed."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GWAS Catalog association file not found: {path}")
    if path.suffix == ".zip":
        archive = zipfile.ZipFile(path)
        members = archive.namelist()
        if len(members) != 1:
            raise ValueError(
                f"{path} should contain exactly the association TSV; found {members}"
            )
        return io.TextIOWrapper(archive.open(members[0]), encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def parse_gwas_associations(
    path: Path,
) -> Tuple[Dict[str, Set[str]], GWASFilterCounts]:
    """Parse the GWAS Catalog file into ``{gene symbol: {trait CURIE}}``.

    Applies the module's row policy (genome-wide significance, intergenic
    drop, mapped-gene column); the returned counts audit every drop.
    """
    logger.info(f"Parsing GWAS Catalog associations from {path}")
    gene_terms: Dict[str, Set[str]] = defaultdict(set)
    n_rows = n_below = n_intergenic = n_flanking = n_no_trait = n_kept = 0

    with _open_association_table(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(_REQUIRED_COLUMNS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"{path} lacks expected column(s) {sorted(missing)}; "
                f"header was {reader.fieldnames}"
            )
        for row in reader:
            n_rows += 1
            try:
                p_value = float(row["P-VALUE"])
            except (TypeError, ValueError):
                n_below += 1
                continue
            if p_value >= GENOME_WIDE_SIGNIFICANCE:
                n_below += 1
                continue
            if (row["INTERGENIC"] or "").strip() == "1":
                n_intergenic += 1
                continue
            mapped_gene = (row["MAPPED_GENE"] or "").strip()
            if not mapped_gene or " - " in mapped_gene:
                n_flanking += 1
                continue
            trait_uris = (row["MAPPED_TRAIT_URI"] or "").strip()
            if not trait_uris:
                n_no_trait += 1
                continue
            traits = {
                trait_uri_to_curie(uri.strip())
                for uri in trait_uris.split(",")
                if uri.strip()
            }
            genes = {
                gene.strip()
                for gene in mapped_gene.replace(";", ",").split(",")
                if gene.strip()
            }
            if not traits or not genes:
                n_no_trait += 1
                continue
            n_kept += 1
            for gene in genes:
                gene_terms[gene] |= traits

    counts = GWASFilterCounts(
        n_rows=n_rows,
        n_below_significance=n_below,
        n_intergenic=n_intergenic,
        n_flanking=n_flanking,
        n_no_trait=n_no_trait,
        n_kept=n_kept,
    )
    logger.info(
        f"  Rows: {counts.n_rows:,}; below genome-wide significance "
        f"(p >= {GENOME_WIDE_SIGNIFICANCE}): {counts.n_below_significance:,}; "
        f"intergenic: {counts.n_intergenic:,} flagged + {counts.n_flanking:,} "
        f"flanking/empty gene; no mapped trait: {counts.n_no_trait:,}; "
        f"kept: {counts.n_kept:,}"
    )
    logger.info(
        f"  Genes: {len(gene_terms):,}; traits: "
        f"{len({term for terms in gene_terms.values() for term in terms}):,}"
    )
    return dict(gene_terms), counts


class GWASCatalogAnnotationSource(AnnotationSource):
    """Domain annotations keyed by GWAS Catalog trait (EFO term).

    Reads the ontology-annotated association file and re-keys its mapped gene
    symbols to UniProt accessions (``DR HGNC`` symbol index of the Swiss-Prot
    flat file) before the statistics see them.
    """

    spec = EFO_SPEC

    def __init__(self, associations_path: Path, dat_path: Path) -> None:
        self.associations_path = Path(associations_path)
        self.dat_path = Path(dat_path)
        #: Populated by :meth:`parse`: the row-policy audit.
        self.filter_counts: Optional[GWASFilterCounts] = None
        #: Populated by :meth:`parse`; its *values* range over gene symbols,
        #: its *keys* over trait CURIEs (axis-swapped, see
        #: remap_gene_annotations).
        self.coverage: Optional[RemapCoverage] = None

    def parse(self) -> Dict[str, Set[str]]:
        gene_terms, self.filter_counts = parse_gwas_associations(self.associations_path)
        index = parse_gene_accession_index(self.dat_path)
        remapped, self.coverage = remap_gene_annotations(
            gene_terms, index.symbol, label="gene symbol→UniProt (GWAS/EFO)"
        )
        return remapped
