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
    enable_true_path: bool
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
            enable_true_path=args.enable_true_path,
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
        missing_inputs=tuple(
            missing_inputs(entry, paths, for_hierarchy=request.enable_true_path)
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
