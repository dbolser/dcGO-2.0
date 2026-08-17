"""Registry of the ontologies the pipeline can associate domains with.

``run_dcgo_human.py --ontology <key>`` used to grow an ``if/elif`` chain per
ontology: one branch to build the :class:`~src.annotation_source.AnnotationSource`,
another to pick the True Path hierarchy. That does not scale to the breadth
UniProt actually offers, so both live here instead — one :class:`OntologyEntry`
per ontology, and the runner just looks the key up.

**What is reachable, and why these.** UniProt is the protein universe for this
pipeline (``protein2ipr``, GOA and Expasy ENZYME are all keyed by UniProt
accession), so any vocabulary the Swiss-Prot flat file already carries per
accession needs *no* identifier mapping. A survey of the human subset of
``uniprot_sprot.dat.gz`` (see ``docs/uniprot_ontology_survey.md``) sorts its
~150 ``DR`` databases into three groups:

* **Vocabularies** — many proteins share a term (Reactome, keywords, OMIM
  phenotypes, Orphanet, TCDB, MEROPS, CAZy, UniPathway, ComplexPortal,
  DrugBank, Pharos, CD-CODE). These are the ones worth testing for domain
  enrichment, and they are registered below.
* **1:1 accession mirrors** — the "term" is just another id for the same
  protein (AlphaFoldDB, STRING, GeneCards, DisGeNET, PhosphoSitePlus …). A
  Fisher test against those is meaningless, so they are not registered; the
  generic ``xref`` escape hatch can still reach them.
* **Domain databases** — Pfam, PANTHER, SUPFAM, Gene3D, CDD, PROSITE, InterPro
  itself. Associating domains with domains is circular, so they are excluded by
  design (again, ``xref`` remains available for deliberate control experiments).

Three further layers are curated into the entry *body* rather than ``DR`` lines
and so need their own extractors (in
:mod:`src.uniprot_annotation_source`): subcellular location, ChEBI ligands /
cofactors, and Rhea reactions.

A fourth kind of entry re-keys one of those vocabularies into a *different term
space* before the statistics see it. ``doid`` / ``orphanet_doid`` take UniProt's
OMIM and Orphanet disease cross-references and translate them to Disease
Ontology terms (:mod:`src.disease_ontology`), which is what gives the disease
layer a hierarchy at all. Note this is a mapping of *terms*, not of proteins —
the protein key space is still UniProt accessions.

Adding an ontology is now a single :class:`OntologyEntry` — the Fisher/FDR
engine and the True Path machinery are untouched.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional

from src.annotation_source import (
    GO_SPEC,
    AnnotationSource,
    GAFAnnotationSource,
    OntologySpec,
)
from src.disease_ontology import DOID_SPEC, DiseaseOntologyAnnotationSource
from src.ec_annotation_source import EC_SPEC, ECAnnotationSource, ec_ancestors
from src.hierarchy import (
    alpha_prefix_ancestors,
    closure_ancestors,
    dotted_ancestors,
    nearest_parents,
    parents_from_map,
    parse_obo_child_parents,
)
from src.relative_inference import ParentsFn
from src.uniprot_annotation_source import (
    COFACTOR_SPEC,
    DISEASE_SPEC,
    KEYWORD_SPEC,
    LIGAND_SPEC,
    REACTOME_SPEC,
    RHEA_SPEC,
    SUBCELLULAR_SPEC,
    UniProtCofactorAnnotationSource,
    UniProtCrossRefAnnotationSource,
    UniProtKeywordAnnotationSource,
    UniProtLigandAnnotationSource,
    UniProtRheaAnnotationSource,
    UniProtSubcellularAnnotationSource,
    parse_keyword_hierarchy,
    parse_reactome_relations,
    parse_subcell_vocabulary,
)

#: An ancestors function: term id → its ancestor term ids.
AncestorsFn = Callable[[str], Iterable[str]]


@dataclass(frozen=True)
class OntologyEntry:
    """One selectable ontology: how to build its annotations and its hierarchy.

    Attributes:
        key: the ``--ontology`` value.
        spec: ontology metadata (id, name, term prefix).
        description: one-line help text.
        build_source: ``(paths, options) → AnnotationSource``.
        build_ancestors: ``(paths) → AncestorsFn``, or ``None`` when no
            hierarchy is available (the run then has no True Path propagation).
        build_parents: ``(paths) → ParentsFn`` giving a term's *direct* parents,
            which is what relative inference ranges over (ancestors are the
            transitive closure and are not a substitute). ``None`` when the
            ontology has no hierarchy, so relative inference is unavailable.
        needs: keys of ``paths`` the source requires, checked before the run so
            a missing input fails loudly instead of half-way through.
        hierarchy_needs: extra ``paths`` keys the hierarchy requires.
        external_propagation: the hierarchy is not handled here. Only GO sets
            this: it propagates through
            :class:`~src.ontology_processor.OntologyProcessor`, which also does
            parental-background (optimal-level) filtering.
    """

    key: str
    spec: OntologySpec
    description: str
    build_source: Callable[[Dict[str, Path], Dict[str, object]], AnnotationSource]
    build_ancestors: Optional[Callable[[Dict[str, Path]], AncestorsFn]] = None
    build_parents: Optional[Callable[[Dict[str, Path]], ParentsFn]] = None
    needs: tuple = ()
    hierarchy_needs: tuple = ()
    external_propagation: bool = False

    @property
    def supports_true_path(self) -> bool:
        return self.build_ancestors is not None or self.external_propagation

    @property
    def supports_relative_inference(self) -> bool:
        """Relative inference needs direct parents, not just ancestors."""
        return self.build_parents is not None or self.external_propagation


_CHILD_PARENTS_CACHE: Dict[tuple, Dict[str, set]] = {}


def _child_parents(
    name: str, path: Path, loader: Callable[[Path], Dict[str, set]]
) -> Dict[str, set]:
    """Parse *path* into a child→parents map once per process.

    Both ``build_ancestors`` and ``build_parents`` derive from the same map, and
    for ChEBI that parse is expensive enough that doing it twice per run is
    worth avoiding. Keyed by (loader name, path) so distinct files never share
    an entry.
    """
    key = (name, str(path))
    if key not in _CHILD_PARENTS_CACHE:
        _CHILD_PARENTS_CACHE[key] = loader(path)
    return _CHILD_PARENTS_CACHE[key]


def _reactome_child_parents(paths: Dict[str, Path]) -> Dict[str, set]:
    return _child_parents(
        "reactome", paths["reactome_relations"], parse_reactome_relations
    )


def _keyword_child_parents(paths: Dict[str, Path]) -> Dict[str, set]:
    return _child_parents("keyword", paths["keywlist"], parse_keyword_hierarchy)


def _dr(
    database: str,
    spec: OntologySpec,
    *,
    id_type: Optional[str] = None,
    term_from_id_type: bool = False,
) -> Callable[[Dict[str, Path], Dict[str, object]], AnnotationSource]:
    """Source factory for a UniProt ``DR`` cross-reference vocabulary."""

    def build(paths: Dict[str, Path], options: Dict[str, object]) -> AnnotationSource:
        return UniProtCrossRefAnnotationSource(
            paths["uniprot_dat"],
            database,
            spec,
            id_type=id_type,
            term_from_id_type=term_from_id_type,
        )

    return build


def _chebi_child_parents(paths: Dict[str, Path]) -> Dict[str, set]:
    return _child_parents("chebi", paths["chebi_obo"], parse_obo_child_parents)


def _chebi_ancestors(paths: Dict[str, Path]) -> AncestorsFn:
    """ChEBI ``is_a`` closure, for the ligand and cofactor layers."""
    return closure_ancestors(_chebi_child_parents(paths))


def _chebi_parents(paths: Dict[str, Path]) -> ParentsFn:
    return parents_from_map(_chebi_child_parents(paths))


def _doid_ancestors(paths: Dict[str, Path]) -> AncestorsFn:
    """Disease Ontology ``is_a`` closure, for the DOID-keyed disease layers.

    ``doid.obo`` carries no ``relationship:`` lines at all — the disease
    classification is pure ``is_a`` — so no extra relations are traversed.
    Obsolete stanzas are excluded (the default), which is what makes the
    ``replaced_by`` resolution in :mod:`src.disease_ontology` necessary.
    """
    return closure_ancestors(_doid_child_parents(paths))


def _doid_child_parents(paths: Dict[str, Path]) -> Dict[str, set]:
    return _child_parents(
        "doid",
        paths["doid_obo"],
        lambda path: parse_obo_child_parents(path, relations=()),
    )


def _doid_parents(paths: Dict[str, Path]) -> ParentsFn:
    return parents_from_map(_doid_child_parents(paths))


def _subcell_child_parents(paths: Dict[str, Path]) -> Dict[str, set]:
    return _child_parents(
        "subcell",
        paths["subcell"],
        lambda path: parse_subcell_vocabulary(path).child_to_parents,
    )


def _subcell_ancestors(paths: Dict[str, Path]) -> AncestorsFn:
    """Subcellular-location closure over ``subcell.txt`` HI/HP edges."""
    return closure_ancestors(_subcell_child_parents(paths))


def _subcell_parents(paths: Dict[str, Path]) -> ParentsFn:
    return parents_from_map(_subcell_child_parents(paths))


def _build_go_source(
    paths: Dict[str, Path], options: Dict[str, object]
) -> AnnotationSource:
    return GAFAnnotationSource(
        paths["gaf"],
        evidence_filter=str(options.get("evidence_filter", "manual")),
        aspects={"P", "F", "C"},
    )


def _build_xref_source(
    paths: Dict[str, Path], options: Dict[str, object]
) -> AnnotationSource:
    """Generic escape hatch: any DR database named at the command line."""
    database = str(options["xref_db"])
    xref_type = options.get("xref_type")
    return UniProtCrossRefAnnotationSource(
        paths["uniprot_dat"],
        database,
        OntologySpec(ontology_id=database, name=f"UniProt {database} cross-reference"),
        id_type=str(xref_type) if xref_type else None,
        term_from_id_type=bool(options.get("xref_term_from_type", False)),
    )


ONTOLOGIES: Dict[str, OntologyEntry] = {
    # ---- Gene Ontology (the reference path) -------------------------------
    "go": OntologyEntry(
        key="go",
        spec=GO_SPEC,
        description="Gene Ontology, from the species GOA GAF file",
        build_source=_build_go_source,
        needs=("gaf",),
        hierarchy_needs=("go_obo",),
        external_propagation=True,
    ),
    # ---- Enzyme Commission (Expasy ENZYME, UniProt-keyed) -----------------
    "ec": OntologyEntry(
        key="ec",
        spec=EC_SPEC,
        description="Enzyme Commission numbers, from Expasy enzyme.dat",
        build_source=lambda paths, options: ECAnnotationSource(paths["enzyme_dat"]),
        build_ancestors=lambda paths: ec_ancestors,
        build_parents=lambda paths: nearest_parents(ec_ancestors),
        needs=("enzyme_dat",),
    ),
    # ---- UniProt DR vocabularies ------------------------------------------
    "reactome": OntologyEntry(
        key="reactome",
        spec=REACTOME_SPEC,
        description="Reactome pathways (DR Reactome)",
        build_source=_dr("Reactome", REACTOME_SPEC),
        build_ancestors=lambda paths: closure_ancestors(_reactome_child_parents(paths)),
        build_parents=lambda paths: parents_from_map(_reactome_child_parents(paths)),
        needs=("uniprot_dat",),
        hierarchy_needs=("reactome_relations",),
    ),
    "keyword": OntologyEntry(
        key="keyword",
        spec=KEYWORD_SPEC,
        description="UniProt keywords (KW lines)",
        build_source=lambda paths, options: UniProtKeywordAnnotationSource(
            paths["uniprot_dat"]
        ),
        build_ancestors=lambda paths: closure_ancestors(_keyword_child_parents(paths)),
        build_parents=lambda paths: parents_from_map(_keyword_child_parents(paths)),
        needs=("uniprot_dat",),
        hierarchy_needs=("keywlist",),
    ),
    "disease": OntologyEntry(
        key="disease",
        spec=DISEASE_SPEC,
        description="OMIM disease phenotypes (DR MIM, phenotype-typed); flat",
        build_source=_dr("MIM", DISEASE_SPEC, id_type="phenotype"),
        needs=("uniprot_dat",),
    ),
    # The same UniProt disease curation, re-keyed onto Disease Ontology terms at
    # parse time. Kept as its own key rather than replacing 'disease': the two
    # test different hypothesis universes (DO pools OMIM's per-locus entries and
    # drops what it does not cross-reference), so comparing them requires both
    # to stay runnable. See src/disease_ontology.py for the mapping policy.
    "doid": OntologyEntry(
        key="doid",
        spec=DOID_SPEC,
        description="Disease Ontology terms, re-keyed from DR MIM at parse time",
        build_source=lambda paths, options: DiseaseOntologyAnnotationSource(
            paths["uniprot_dat"], paths["doid_obo"]
        ),
        build_ancestors=_doid_ancestors,
        build_parents=_doid_parents,
        needs=("uniprot_dat", "doid_obo"),
        hierarchy_needs=("doid_obo",),
    ),
    "orphanet": OntologyEntry(
        key="orphanet",
        spec=OntologySpec(ontology_id="Orphanet", name="Orphanet rare disease"),
        description="Orphanet rare diseases (DR Orphanet); flat",
        build_source=_dr(
            "Orphanet",
            OntologySpec(ontology_id="Orphanet", name="Orphanet rare disease"),
        ),
        needs=("uniprot_dat",),
    ),
    # DO cross-references Orphanet too (xref: ORDO:<id>), so the same machinery
    # gives the Orphanet layer a hierarchy for the price of a different prefix.
    "orphanet_doid": OntologyEntry(
        key="orphanet_doid",
        spec=DOID_SPEC,
        description="Disease Ontology terms, re-keyed from DR Orphanet at parse time",
        build_source=lambda paths, options: DiseaseOntologyAnnotationSource(
            paths["uniprot_dat"],
            paths["doid_obo"],
            database="Orphanet",
            id_type=None,
            xref_prefix="ORDO",
        ),
        build_ancestors=_doid_ancestors,
        build_parents=_doid_parents,
        needs=("uniprot_dat", "doid_obo"),
        hierarchy_needs=("doid_obo",),
    ),
    "tcdb": OntologyEntry(
        key="tcdb",
        spec=OntologySpec(
            ontology_id="TCDB", name="Transporter Classification Database"
        ),
        description="Transporter classification, e.g. 8.A.98.1.10 (DR TCDB)",
        build_source=_dr(
            "TCDB",
            OntologySpec(
                ontology_id="TCDB", name="Transporter Classification Database"
            ),
        ),
        # TC numbers nest exactly like EC: class.subclass.family.subfamily.system.
        build_ancestors=lambda paths: dotted_ancestors,
        build_parents=lambda paths: nearest_parents(dotted_ancestors),
        needs=("uniprot_dat",),
    ),
    "merops": OntologyEntry(
        key="merops",
        spec=OntologySpec(ontology_id="MEROPS", name="MEROPS peptidase classification"),
        description="Peptidase/inhibitor families, e.g. S01.151 (DR MEROPS)",
        build_source=_dr(
            "MEROPS",
            OntologySpec(ontology_id="MEROPS", name="MEROPS peptidase classification"),
        ),
        # "S01.151" → family "S01" → catalytic type "S".
        build_ancestors=lambda paths: alpha_prefix_ancestors,
        build_parents=lambda paths: nearest_parents(alpha_prefix_ancestors),
        needs=("uniprot_dat",),
    ),
    "cazy": OntologyEntry(
        key="cazy",
        spec=OntologySpec(ontology_id="CAZy", name="Carbohydrate-Active enZymes"),
        description="CAZy families, e.g. GT32 (DR CAZy)",
        build_source=_dr(
            "CAZy", OntologySpec(ontology_id="CAZy", name="Carbohydrate-Active enZymes")
        ),
        # "GT32" → class "GT" (glycosyltransferase).
        build_ancestors=lambda paths: alpha_prefix_ancestors,
        build_parents=lambda paths: nearest_parents(alpha_prefix_ancestors),
        needs=("uniprot_dat",),
    ),
    "unipathway": OntologyEntry(
        key="unipathway",
        spec=OntologySpec(
            ontology_id="UniPathway", name="UniPathway", term_prefix="UPA"
        ),
        description="UniPathway metabolic pathways (DR UniPathway)",
        build_source=_dr(
            "UniPathway",
            OntologySpec(
                ontology_id="UniPathway", name="UniPathway", term_prefix="UPA"
            ),
        ),
        needs=("uniprot_dat",),
    ),
    "complex": OntologyEntry(
        key="complex",
        spec=OntologySpec(
            ontology_id="ComplexPortal", name="Protein complex", term_prefix="CPX-"
        ),
        description="Protein complexes (DR ComplexPortal)",
        build_source=_dr(
            "ComplexPortal",
            OntologySpec(
                ontology_id="ComplexPortal", name="Protein complex", term_prefix="CPX-"
            ),
        ),
        needs=("uniprot_dat",),
    ),
    "drugbank": OntologyEntry(
        key="drugbank",
        spec=OntologySpec(
            ontology_id="DrugBank", name="DrugBank drug", term_prefix="DB"
        ),
        description="DrugBank drugs targeting the protein (DR DrugBank)",
        build_source=_dr(
            "DrugBank",
            OntologySpec(
                ontology_id="DrugBank", name="DrugBank drug", term_prefix="DB"
            ),
        ),
        needs=("uniprot_dat",),
    ),
    "pharos": OntologyEntry(
        key="pharos",
        spec=OntologySpec(ontology_id="Pharos", name="Pharos target development level"),
        # DR Pharos; P31946; Tbio.  — the id is the accession, the vocabulary
        # (Tclin/Tchem/Tbio/Tdark) is the third field.
        description="Pharos target development level: Tclin/Tchem/Tbio/Tdark",
        build_source=_dr(
            "Pharos",
            OntologySpec(ontology_id="Pharos", name="Pharos target development level"),
            term_from_id_type=True,
        ),
        needs=("uniprot_dat",),
    ),
    "condensate": OntologyEntry(
        key="condensate",
        spec=OntologySpec(ontology_id="CD-CODE", name="Biomolecular condensate"),
        # DR CD-CODE; 91857CE7; Nucleolus.  — again the name is the third field.
        description="Biomolecular condensates, e.g. Nucleolus (DR CD-CODE)",
        build_source=_dr(
            "CD-CODE",
            OntologySpec(ontology_id="CD-CODE", name="Biomolecular condensate"),
            term_from_id_type=True,
        ),
        needs=("uniprot_dat",),
    ),
    # ---- UniProt entry-body layers (CC comments and FT qualifiers) --------
    "subcellular": OntologyEntry(
        key="subcellular",
        spec=SUBCELLULAR_SPEC,
        description="Subcellular locations (CC SUBCELLULAR LOCATION → SL- terms)",
        build_source=lambda paths, options: UniProtSubcellularAnnotationSource(
            paths["uniprot_dat"], paths["subcell"]
        ),
        build_ancestors=_subcell_ancestors,
        build_parents=_subcell_parents,
        needs=("uniprot_dat", "subcell"),
        hierarchy_needs=("subcell",),
    ),
    "ligand": OntologyEntry(
        key="ligand",
        spec=LIGAND_SPEC,
        description="Bound ligands from FT /ligand_id ChEBI qualifiers",
        build_source=lambda paths, options: UniProtLigandAnnotationSource(
            paths["uniprot_dat"]
        ),
        build_ancestors=_chebi_ancestors,
        build_parents=_chebi_parents,
        needs=("uniprot_dat",),
        hierarchy_needs=("chebi_obo",),
    ),
    "cofactor": OntologyEntry(
        key="cofactor",
        spec=COFACTOR_SPEC,
        description="Cofactors from CC COFACTOR ChEBI cross-references",
        build_source=lambda paths, options: UniProtCofactorAnnotationSource(
            paths["uniprot_dat"]
        ),
        build_ancestors=_chebi_ancestors,
        build_parents=_chebi_parents,
        needs=("uniprot_dat",),
        hierarchy_needs=("chebi_obo",),
    ),
    "rhea": OntologyEntry(
        key="rhea",
        spec=RHEA_SPEC,
        description="Catalysed reactions from CC CATALYTIC ACTIVITY (Rhea)",
        build_source=lambda paths, options: UniProtRheaAnnotationSource(
            paths["uniprot_dat"]
        ),
        needs=("uniprot_dat",),
    ),
    # ---- Escape hatch ------------------------------------------------------
    "xref": OntologyEntry(
        key="xref",
        spec=OntologySpec(ontology_id="xref", name="Arbitrary UniProt DR database"),
        description="Any DR database named by --xref-db (e.g. KEGG, BRENDA)",
        build_source=_build_xref_source,
        needs=("uniprot_dat",),
    ),
}


def ontology_keys() -> List[str]:
    """Selectable ``--ontology`` values, in registration order."""
    return list(ONTOLOGIES)


def get_ontology(key: str) -> OntologyEntry:
    """Look up an ontology by ``--ontology`` key.

    Raises:
        KeyError: with the list of valid keys, for an unknown ontology.
    """
    try:
        return ONTOLOGIES[key]
    except KeyError:
        raise KeyError(
            f"Unknown ontology {key!r}. Available: {', '.join(ontology_keys())}"
        ) from None


def describe_ontologies() -> str:
    """Multi-line ``key — description`` listing, for CLI help."""
    width = max(len(key) for key in ONTOLOGIES)
    return "\n".join(
        f"  {entry.key:<{width}}  {entry.description}" for entry in ONTOLOGIES.values()
    )


def missing_inputs(
    entry: OntologyEntry, paths: Dict[str, Path], *, for_hierarchy: bool = False
) -> List[str]:
    """Return the ``paths`` entries this ontology needs but that are absent.

    ``for_hierarchy`` additionally checks the True Path inputs. An empty list
    means the run can proceed.
    """
    required = entry.needs + (entry.hierarchy_needs if for_hierarchy else ())
    missing = []
    for name in required:
        path = paths.get(name)
        if path is None or not Path(path).exists():
            missing.append(f"{name} ({path})" if path is not None else name)
    return missing
