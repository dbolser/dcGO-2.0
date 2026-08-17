"""Species- and ontology-generic public runner boundary.

The root :mod:`run_dcgo_human` module remains the implementation and script
compatibility entry point while orchestration stages move here incrementally.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from src.ontology_registry import OntologyEntry, get_ontology, missing_inputs


@dataclass(frozen=True, slots=True)
class RunRequest:
    """Typed, growing subset of options consumed by extracted run stages."""

    species: str
    ontology: str
    domain_key: str
    evidence_filter: str
    output_dir: Path
    fdr_threshold: float
    min_support: int
    min_ic: float
    enable_true_path: bool
    enable_relative_inference: bool
    propagate_annotations: bool
    enable_supra_domains: bool
    xref_db: str | None
    go_ontology: Path
    enzyme_dat: Path
    uniprot_dat: Path
    reactome_relations: Path
    keyword_list: Path
    subcell: Path
    chebi_obo: Path
    doid_obo: Path
    hpo_g2p: Path
    hpo_obo: Path
    syngo_zip: Path
    mgi_genepheno: Path
    mgi_marker_swissprot: Path
    mp_obo: Path
    wb_phenotype: Path
    worm_idmapping: Path
    wbphenotype_obo: Path
    zfin_phenotype: Path
    zfin_uniprot: Path
    zfa_obo: Path
    fb_genotype_phenotype: Path
    fbal_to_fbgn: Path
    fbgn_uniprot: Path
    fbbt_obo: Path
    fbcv_obo: Path

    @property
    def engages_hierarchy(self) -> bool:
        """Whether any stage of this run reads the ontology's hierarchy.

        True for input-map propagation (``--propagate-annotations``), either
        Stage 5.5/4.5 hierarchy stage (``--enable-true-path``,
        ``--enable-relative-inference``), and the ``--min-ic`` floor, whose
        frequencies are estimated over the propagated map. The single
        authority for "this run needs the hierarchy inputs": input
        resolution, the manifest's hierarchy_inputs record, and the ic_source
        decision all read this predicate.
        """
        return (
            self.propagate_annotations
            or self.enable_true_path
            or self.enable_relative_inference
            or self.min_ic > 0
        )

    @classmethod
    def from_namespace(cls, args: argparse.Namespace) -> "RunRequest":
        """Translate the compatibility parser's namespace into generic names."""
        return cls(
            species=args.species,
            ontology=args.ontology,
            domain_key=getattr(args, "domain_key", "interpro"),
            evidence_filter=args.evidence_filter,
            output_dir=args.output_dir,
            fdr_threshold=args.fdr_threshold,
            min_support=args.min_support,
            min_ic=getattr(args, "min_ic", 0.0),
            enable_true_path=args.enable_true_path,
            enable_relative_inference=getattr(args, "enable_relative_inference", False),
            propagate_annotations=getattr(args, "propagate_annotations", False),
            enable_supra_domains=args.enable_supra_domains,
            xref_db=getattr(args, "xref_db", None),
            go_ontology=args.go_ontology,
            enzyme_dat=args.enzyme_dat,
            uniprot_dat=args.uniprot_dat,
            reactome_relations=args.reactome_relations,
            keyword_list=args.keyword_list,
            subcell=args.subcell,
            chebi_obo=args.chebi_obo,
            doid_obo=args.doid_obo,
            hpo_g2p=args.hpo_genes_to_phenotype,
            hpo_obo=args.hpo_obo,
            syngo_zip=args.syngo_zip,
            mgi_genepheno=args.mgi_genepheno,
            mgi_marker_swissprot=args.mgi_marker_swissprot,
            mp_obo=args.mp_obo,
            wb_phenotype=args.wb_phenotype,
            worm_idmapping=args.worm_idmapping,
            wbphenotype_obo=args.wbphenotype_obo,
            zfin_phenotype=args.zfin_phenotype,
            zfin_uniprot=args.zfin_uniprot,
            zfa_obo=args.zfa_obo,
            fb_genotype_phenotype=args.fb_genotype_phenotype,
            fbal_to_fbgn=args.fbal_to_fbgn,
            fbgn_uniprot=args.fbgn_uniprot,
            fbbt_obo=args.fbbt_obo,
            fbcv_obo=args.fbcv_obo,
        )

    def ontology_paths(self) -> dict[str, Path]:
        """Return every registry-addressable ontology input path."""
        return {
            "gaf": Path(f"data/raw/goa_annotations/goa_{self.species}.gaf.gz"),
            "go_obo": self.go_ontology,
            "enzyme_dat": self.enzyme_dat,
            "uniprot_dat": self.uniprot_dat,
            "reactome_relations": self.reactome_relations,
            "keywlist": self.keyword_list,
            "subcell": self.subcell,
            "chebi_obo": self.chebi_obo,
            "doid_obo": self.doid_obo,
            "hpo_g2p": self.hpo_g2p,
            "hpo_obo": self.hpo_obo,
            "syngo_zip": self.syngo_zip,
            "mgi_genepheno": self.mgi_genepheno,
            "mgi_marker_swissprot": self.mgi_marker_swissprot,
            "mp_obo": self.mp_obo,
            "wb_phenotype": self.wb_phenotype,
            "worm_idmapping": self.worm_idmapping,
            "wbphenotype_obo": self.wbphenotype_obo,
            "zfin_phenotype": self.zfin_phenotype,
            "zfin_uniprot": self.zfin_uniprot,
            "zfa_obo": self.zfa_obo,
            "fb_genotype_phenotype": self.fb_genotype_phenotype,
            "fbal_to_fbgn": self.fbal_to_fbgn,
            "fbgn_uniprot": self.fbgn_uniprot,
            "fbbt_obo": self.fbbt_obo,
            "fbcv_obo": self.fbcv_obo,
        }


@dataclass(frozen=True, slots=True)
class InputResolution:
    """Registry identity and validated paths needed before expensive stages."""

    ontology_entry: OntologyEntry
    ontology_label: str
    ontology_paths: dict[str, Path]
    true_path_unsupported: bool
    missing_inputs: tuple[str, ...]


def resolve_inputs(request: RunRequest) -> InputResolution:
    """Resolve the ontology registry entry and all inputs for *request*."""
    entry = get_ontology(request.ontology)
    paths = request.ontology_paths()
    return InputResolution(
        ontology_entry=entry,
        ontology_label=(
            request.xref_db.lower()
            if request.ontology == "xref" and request.xref_db
            else request.ontology
        ),
        ontology_paths=paths,
        true_path_unsupported=(
            request.enable_true_path and not entry.supports_true_path
        ),
        # Anything that engages the hierarchy (see RunRequest.engages_hierarchy)
        # makes the hierarchy inputs mandatory up front.
        missing_inputs=tuple(
            missing_inputs(
                entry,
                paths,
                for_hierarchy=request.engages_hierarchy,
            )
        ),
    )


def parse_run_request(argv: list[str] | None = None) -> RunRequest:
    """Parse *argv* into the stable package-facing request type."""
    from run_dcgo_human import parse_arguments

    args, _parser = parse_arguments(argv)
    return RunRequest.from_namespace(args)


def main(argv: list[str] | None = None) -> int:
    """Run domain/ontology inference for any supported species and ontology."""
    from run_dcgo_human import main as legacy_main

    return legacy_main(argv)
