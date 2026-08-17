"""ZFIN zebrafish phenotypes as a ZFA (anatomy) layer, re-keyed to UniProt.

The dcGO trick applied to zebrafish (``--species zebrafish --ontology zfa``):
domain → ZFA associations are learned on zebrafish proteins, and the
species-agnostic domains carry them anywhere.

**Semantics — read this before interpreting the layer.** ZFIN curates
phenotypes as EQ post-compositions (Entity = an anatomical structure or
process, Quality = a PATO term), not as pre-composed phenotype-ontology ids.
This layer therefore does *not* say "domain D is associated with phenotype P";
it says **"proteins whose mutation produces an abnormality of anatomical
structure S are enriched for domain D"** — gene → anatomical-structure-affected,
propagated up the ZFA ``is_a``/``part_of`` DAG.

Why not ZP, the pre-composed Zebrafish Phenotype ontology? Assessed and
rejected as not derivable from what is on disk: ``phenoGeneCleanData_fish.txt``
carries no ZP ids, and composing EQ → ZP needs each ZP term's logical
definition, which the OBO edition of ``zp.obo`` does not carry (0 of 43,521
non-obsolete ZP terms have ``intersection_of`` lines; only 924 carry a
machine-readable EQ signature in a ``comment:``, i.e. ~2% — the definitions
live only in the OWL equivalence axioms). Falling back to the affected-anatomy
entity is honest and loses only the quality dimension: 169,887 of 169,887 rows
at acquisition are tagged ``abnormal``, so "some abnormality of S" is exactly
what the rows assert.

Parsing policy for ``phenoGeneCleanData_fish.txt`` (25 unnamed tab columns):
the gene is column 3 (``ZDB-GENE-``); the affected entity is taken from *all
four* entity columns — E1 subterm (4), E1 superterm (8), E2 subterm (13), E2
superterm (17) — because sub- and superterm are both genuinely affected
structures in a post-composition ("cell X of tissue Y"). Only ``ZFA:`` ids are
kept: the entity slots also carry GO process ids (~46k at acquisition), BSPO
positional and CHEBI ids, which do not belong to an anatomy layer and are
counted as skipped. Rows not tagged ``abnormal`` would be dropped and counted
(none exist in the current file). The clean-data file already excludes
morphants and multi-gene fish — that filtering is ZFIN's, not ours.

The ZDB-GENE → UniProt translation is ZFIN's own ``uniprot.txt`` (gene id,
SO type, symbol, accession — one row per pairing, TrEMBL included), applied at
parse time with the counted unmapped/one-to-many policy
(:func:`src.gene_mapping.remap_gene_annotations`).

The hierarchy is ``zfa.obo``, read by the shared light OBO reader.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Set

from loguru import logger

from src.annotation_source import AnnotationSource, OntologySpec
from src.gene_mapping import GeneAccessionMap, remap_gene_annotations
from src.hierarchy import open_text
from src.remap import RemapCoverage

#: Zebrafish Anatomy Ontology terms, e.g. ``ZFA:0000108``.
ZFA_SPEC = OntologySpec(
    ontology_id="ZFA", name="Zebrafish Anatomy Ontology", term_prefix="ZFA:"
)

#: Column indices of phenoGeneCleanData_fish.txt (no header row).
_GENE_COL = 2
#: E1 subterm, E1 superterm, E2 subterm, E2 superterm.
_ENTITY_COLS = (3, 7, 12, 16)
_TAG_COL = 11


def parse_pheno_gene_clean_data(path: Path) -> Dict[str, Set[str]]:
    """Parse ``phenoGeneCleanData_fish.txt`` into ``{ZDB-GENE id: {ZFA term}}``.

    See the module docstring for which columns are read and why non-``ZFA:``
    entities are skipped.
    """
    logger.info(f"Parsing ZFIN phenotype annotations from {path}")
    gene_terms: Dict[str, Set[str]] = defaultdict(set)
    n_rows = 0
    n_not_abnormal = 0
    n_non_gene = 0
    n_entities = 0
    n_non_zfa = 0
    with open_text(path, label="ZFIN phenotype file") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(_ENTITY_COLS):
                continue
            gene = fields[_GENE_COL].strip()
            if not gene.startswith("ZDB-GENE-"):
                # ZFIN also phenotypes non-gene loci (ncRNA classes such as
                # ZDB-MIRNAG-, ZDB-LINCRNAG-). They have no protein product to
                # re-key, so they are dropped — counted like every other drop.
                n_non_gene += 1
                continue
            if fields[_TAG_COL].strip() != "abnormal":
                n_not_abnormal += 1
                continue
            n_rows += 1
            for col in _ENTITY_COLS:
                entity = fields[col].strip()
                if not entity:
                    continue
                if entity.startswith("ZFA:"):
                    n_entities += 1
                    gene_terms[gene].add(entity)
                else:
                    n_non_zfa += 1
    logger.info(
        f"  Abnormal rows: {n_rows:,} (non-abnormal dropped: {n_not_abnormal:,}; "
        f"non-gene loci dropped: {n_non_gene:,}); "
        f"ZFA entities kept: {n_entities:,}; non-ZFA entities skipped "
        f"(GO/BSPO/CHEBI/...): {n_non_zfa:,}"
    )
    logger.info(
        f"  Genes: {len(gene_terms):,}; terms: "
        f"{len({term for terms in gene_terms.values() for term in terms}):,}"
    )
    return dict(gene_terms)


def parse_zfin_uniprot(path: Path) -> GeneAccessionMap:
    """Parse ZFIN's ``uniprot.txt`` into a ZDB-GENE → accessions map.

    The file also cross-references non-gene ZDB objects (BACs etc.); only
    ``ZDB-GENE-`` rows are kept. One gene may list several accessions — kept
    one-to-many for the counted expansion policy downstream.
    """
    logger.info(f"Building ZDB-GENE → accession map from {path}")
    mapping: Dict[str, Set[str]] = defaultdict(set)
    n_rows = 0
    with open_text(path, label="ZFIN uniprot.txt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            gene, accession = fields[0].strip(), fields[3].strip()
            if gene.startswith("ZDB-GENE-") and accession:
                n_rows += 1
                mapping[gene].add(accession)
    gene_map = GeneAccessionMap("ZDB-GENE", dict(mapping), n_rows)
    logger.info(
        f"  Genes with accessions: {len(gene_map):,} "
        f"({gene_map.n_one_to_many:,} one-to-many)"
    )
    return gene_map


class ZFINAnatomyAnnotationSource(AnnotationSource):
    """Domain annotations keyed by affected zebrafish anatomical structure.

    Reads ``phenoGeneCleanData_fish.txt`` (abnormal rows, ZFA entities) and
    re-keys its ZDB-GENE ids to UniProt accessions via ZFIN's ``uniprot.txt``.
    """

    spec = ZFA_SPEC

    def __init__(self, phenotype_path: Path, uniprot_map_path: Path) -> None:
        self.phenotype_path = Path(phenotype_path)
        self.uniprot_map_path = Path(uniprot_map_path)
        #: Populated by :meth:`parse`; its *values* range over ZDB-GENE ids,
        #: its *keys* over ZFA terms (axis-swapped, see remap_gene_annotations).
        self.coverage: Optional[RemapCoverage] = None

    def parse(self) -> Dict[str, Set[str]]:
        gene_terms = parse_pheno_gene_clean_data(self.phenotype_path)
        gene_map = parse_zfin_uniprot(self.uniprot_map_path)
        remapped, self.coverage = remap_gene_annotations(
            gene_terms,
            gene_map,
            label="ZDB-GENE→UniProt (ZFA)",
            target_label="UniProt accession",
        )
        return remapped
