#!/usr/bin/env python3
"""
dcGO Pipeline - Protein Domain-GO Association Analysis

This script runs the complete dcGO statistical inference pipeline for any species.
It performs domain-GO association analysis using sparse matrix operations and parallel
Fisher's exact tests.

Every run writes a provenance manifest, run_manifest_<ontology>.json, into the
output directory: input/output SHA-256 hashes and release headers, the Git
revision, the uv.lock hash, the command line and every effective threshold. See
REPRODUCIBILITY.md.

Usage:
    uv run python run_dcgo_human.py [OPTIONS]

Options:
    --species STR            Species to analyze: 'human', 'mouse', etc. (default: human)
    --ontology STR           Ontology to associate domains with (default: go). See
                             src/ontology_registry.py, or --help, for the full list:
                             go, ec, reactome, keyword, disease, doid, orphanet,
                             orphanet_doid, tcdb, merops, cazy, unipathway, complex,
                             drugbank, pharos, condensate, subcellular, ligand,
                             cofactor, rhea, xref
    --doid-obo PATH          Path to doid.obo, used when --ontology doid|orphanet_doid
    --xref-db STR            UniProt DR database name, required when --ontology xref (e.g. KEGG, BRENDA)
    --xref-type STR          Optional DR third-field filter for --ontology xref (e.g. 'phenotype')
    --enzyme-dat PATH        Path to Expasy enzyme.dat, used when --ontology ec
    --uniprot-dat PATH       Path to UniProt Swiss-Prot flat file, used by every UniProt-native ontology
    --subcell PATH           Path to UniProt subcell.txt, used when --ontology subcellular
    --chebi-obo PATH         Path to ChEBI OBO, for --ontology ligand|cofactor --enable-true-path
    --domain-key STR         Which protein2ipr column defines a domain: 'interpro'
                             (integrated entry, default) or 'ssf' (SUPERFAMILY/SCOP
                             signature, the published dcGO's domain universe)
    --evidence-filter STR    Evidence code filter: 'all', 'manual', 'experimental' (default: manual)
    --fdr-threshold FLOAT    FDR significance threshold (default: 0.01)
    --num-cores INT          Accepted for compatibility; the Fisher stage is in-process Cython (default: 8)
    --output-dir PATH        Output directory for results (default: results/)
    --batch-size INT         Batch size for Fisher tests (default: 50000)
    --permute-annotations N  Calibration control: shuffle protein↔term-set assignment (null run)
    --enable-true-path       True Path propagation (GO via OBO DAG, EC via numbering, reactome/keyword via hierarchy files)
    --go-ontology PATH       Path to GO ontology file (default: data/raw/go_ontology/go-basic.obo)

Examples:
    # Run for human proteins
    uv run python run_dcgo_human.py --species human --num-cores 16

    # Run for mouse proteins with experimental evidence only
    uv run python run_dcgo_human.py --species mouse --evidence-filter experimental

    # Associate human domains with Enzyme Commission numbers instead of GO
    uv run python run_dcgo_human.py --ontology ec

    # UniProt-native vocabularies (all keyed by accession, no id mapping)
    uv run python run_dcgo_human.py --ontology reactome
    uv run python run_dcgo_human.py --ontology keyword
    uv run python run_dcgo_human.py --ontology disease            # OMIM phenotype (DR MIM), flat
    uv run python run_dcgo_human.py --ontology doid --enable-true-path  # same curation, on the DO DAG
    uv run python run_dcgo_human.py --ontology subcellular        # CC SUBCELLULAR LOCATION
    uv run python run_dcgo_human.py --ontology ligand             # FT /ligand_id (ChEBI)
    uv run python run_dcgo_human.py --ontology tcdb               # transporter classification
    uv run python run_dcgo_human.py --ontology xref --xref-db KEGG # any DR database

    # Run with True Path Rule propagation
    uv run python run_dcgo_human.py --enable-true-path --go-ontology data/raw/go_ontology/go-basic.obo

    # Key domains by SCOP superfamily instead of InterPro entry, to compare
    # against the published dcGO (VALIDATION_PLAN §3)
    uv run python run_dcgo_human.py --domain-key ssf
"""

import argparse
import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from loguru import logger
from scipy.stats import hypergeom

from src.annotation_source import restrict_to_universe
from src.domain_annotation_parser import DOMAIN_KEYS, DomainAnnotationParser
from src.hierarchy import propagate_via_ancestors
from src.ontology_processor import OntologyProcessor
from src.ontology_registry import (
    OntologyEntry,
    describe_ontologies,
    get_ontology,
    missing_inputs,
    ontology_keys,
)
from src.run_manifest import RunManifest, describe_file, manifest_filename
from src.sparse_fisher import (
    DomainType,
    build_sparse_matrices,
    compute_contingency_tables_sparse,
)
from src.vectorized_fisher import (
    benjamini_hochberg_by_family,
    fisher_exact_parallel,
)

logger.remove()
logger.add(sys.stderr, level="INFO")

# Analysis constants that are not exposed as CLI flags. They are named here (and
# recorded in the run manifest) rather than buried as literals at their call
# sites, so "every threshold" in the provenance record cannot drift from the
# thresholds the code actually applied.
MAX_SUPRA_DOMAIN_LENGTH = 3  # longest contiguous domain combination tested
MIN_DOMAIN_LENGTH = 10  # residues; shorter InterPro matches are discarded
FISHER_ALTERNATIVE = "greater"  # one-sided: enrichment only
PARENTAL_MIN_BACKGROUND_SIZE = 3  # GO parental-background (optimal-level) test
PARENTAL_ALPHA_THRESHOLD = 0.05  # raw-p cutoff for that parent/child test

#: ``ontology_paths`` key → the ``config/settings.py`` data source it is
#: downloaded from. Used only to *label* manifest inputs with where they came
#: from; the SHA-256 recorded alongside is what actually identifies the bytes.
INPUT_SOURCE_NAMES = {
    "go_obo": "go_ontology",
    "enzyme_dat": "enzyme",
    "uniprot_dat": "uniprot_sprot_dat",
    "reactome_relations": "reactome_relations",
    "keywlist": "uniprot_keywlist",
    "subcell": "uniprot_subcell",
    "chebi_obo": "chebi",
    # Not an ontology input: the upstream source protein2ipr_<species>.dat.gz
    # was derived from, recorded as `derived_from`.
    "interpro_mappings": "interpro_mappings",
}


def input_source_urls(species: str) -> dict:
    """Upstream download URL for each pipeline input, keyed as in ``INPUT_SOURCE_NAMES``.

    Best-effort by design: ``config/`` is a source-checkout convenience and is
    not shipped in the wheel, and importing it builds the project directory
    layout as a side effect. A run that cannot resolve URLs records input hashes
    without them rather than failing — the hash is the identity, the URL is a
    convenience label.
    """
    try:
        from config.settings import Config

        config = Config(use_env_overrides=False)
    except Exception as exc:  # any config problem: these labels are optional
        logger.debug(f"Source URLs unavailable for the manifest: {exc}")
        return {}

    urls = {
        key: config.data_sources[name].url
        for key, name in INPUT_SOURCE_NAMES.items()
        if name in config.data_sources
    }
    # data_sources pins the *human* GAF; every other species has its own URL.
    urls["gaf"] = config.goa_url_for(species)
    return urls


@dataclass
class AssociationResult:
    """Simple dataclass to hold association results for True Path Rule."""

    domain: str
    go_term: str
    p_value: float
    q_value: float
    hyper_score: float
    a: int  # proteins with both
    b: int  # proteins with domain only
    c: int  # proteins with GO only
    d: int  # proteins with neither


def calculate_hypergeometric_score(a: int, b: int, c: int, d: int) -> float:
    """
    Calculate hypergeometric-based association score on 1-100 scale.

    Args:
        a: Proteins with both domain and GO term
        b: Proteins with GO term only (not domain)
        c: Proteins with domain only (not GO)
        d: Proteins with neither

    Returns:
        float: Association score between 1.0 and 100.0, or NaN when the score
        could not be computed. NaN is deliberate: this column is exported and
        read downstream, so a numerical failure must be visibly missing rather
        than a plausible mid-range number. It used to return 50.0 here, which is
        indistinguishable from a genuine medium-confidence association.
    """
    n = a + b + c + d  # total proteins
    k = a + c  # proteins with domain
    m = a + b  # proteins with GO term
    x = a  # proteins with both

    # A contingency cell cannot be negative. Reject rather than hand the values
    # to hypergeom, which returns NaN for them without raising — see below.
    if min(a, b, c, d) < 0:
        logger.warning(
            f"Hypergeometric score: negative contingency cell "
            f"(a={a} b={b} c={c} d={d}). Reporting NaN."
        )
        return float("nan")

    if k == 0 or m == 0 or x == 0:
        return 0.0

    try:
        # Calculate hypergeometric survival function (1 - CDF)
        # P(X ≥ x) where X ~ Hypergeometric(n, k, m)
        p_hyper = hypergeom.sf(x - 1, n, k, m)

        # NaN and zero must not share a branch. `scipy.stats.hypergeom.sf`
        # returns NaN for invalid parameters *without raising*, so folding NaN
        # into the "p is too small to represent" case reported the maximum
        # score, 100.0, for a table that could not be evaluated at all — the
        # loudest possible answer to an unanswerable question.
        if np.isnan(p_hyper):
            logger.warning(
                f"Hypergeometric score is NaN for a={a} b={b} c={c} d={d}. "
                "Reporting NaN."
            )
            return float("nan")

        if p_hyper > 0:
            # Convert to -log10 scale
            score = -np.log10(p_hyper)
            # Scale to 1-100 range (typical values 1e-50 to 1e-1 give scores 1-500)
            return float(min(100.0, max(1.0, score * 10)))

        # A genuine underflow to exactly 0: the association is as strong as this
        # scale can express.
        return 100.0

    except (ValueError, OverflowError, ZeroDivisionError) as exc:
        # Explicitly missing, not "neutral". See the Returns note above.
        logger.warning(
            f"Hypergeometric score failed for a={a} b={b} c={c} d={d}: {exc}. "
            "Reporting NaN."
        )
        return float("nan")


def odds_ratio_interval(a: int, b: int, c: int, d: int) -> tuple:
    """Woolf log-interval for the odds ratio, Haldane–Anscombe corrected.

    FDR significance alone can keep biologically fragile associations built on
    sparse tables — a P1 review item asks for the contingency cells and an
    odds-ratio interval so a reader can see that fragility instead of inferring
    it from a bare p-value.

    Method, stated because it is a choice and not the only one: 0.5 is added to
    every cell when any cell is zero (Haldane–Anscombe), then the interval is
    ``exp(log OR ± 1.96 × sqrt(1/a + 1/b + 1/c + 1/d))``. This is the standard
    large-sample interval; it is *not* the exact conditional interval, so on the
    smallest tables treat it as indicative. `VALIDATION_PLAN.md` §5 carries the
    open decision about adopting the Haldane-corrected OR as the point estimate
    too — this changes no existing column.

    Returns ``(low, high)``, or ``(nan, nan)`` when the table is unusable.
    """
    if min(a, b, c, d) < 0:
        return (float("nan"), float("nan"))
    cells = [a, b, c, d]
    if min(cells) == 0:
        cells = [value + 0.5 for value in cells]
    ca, cb, cc, cd = cells
    if cb == 0 or cc == 0:
        return (float("nan"), float("nan"))
    try:
        log_or = math.log((ca * cd) / (cb * cc))
        se = math.sqrt(1 / ca + 1 / cb + 1 / cc + 1 / cd)
    except (ValueError, ZeroDivisionError, OverflowError):
        return (float("nan"), float("nan"))
    return (math.exp(log_or - 1.96 * se), math.exp(log_or + 1.96 * se))


def validate_arguments(
    args: argparse.Namespace, parser: argparse.ArgumentParser
) -> None:
    """Reject nonsensical parameters before the expensive stages start.

    argparse checks types, not ranges. Without this an ``--fdr-threshold 5``
    past the prior instead of interpolating toward it, and a non-positive
    ``--batch-size`` makes the progress callback divide by zero — each of them an
    hour into a run, or worse, silently. ``parser.error`` exits 2 with usage,
    which is what a caller checking the exit status expects.
    """
    if not 0.0 < args.fdr_threshold <= 1.0:
        parser.error(
            f"--fdr-threshold must be in (0, 1], got {args.fdr_threshold}. "
            "A threshold above 1 calls every test significant."
        )
    if args.batch_size <= 0:
        parser.error(f"--batch-size must be positive, got {args.batch_size}")
    if args.min_support < 0:
        parser.error(f"--min-support must be >= 0, got {args.min_support}")
    if args.num_cores <= 0:
        parser.error(f"--num-cores must be positive, got {args.num_cores}")
    if not args.species or "/" in args.species:
        parser.error(
            f"--species must be a bare name used in the input filenames, got "
            f"{args.species!r}"
        )


def build_ontology_paths(args: argparse.Namespace) -> dict:
    """Resolve every input any ontology might need, keyed as the registry names them.

    One place, so the registry (``src/ontology_registry.py``) is free to say
    which of these a given ontology actually uses for its annotations and for
    its hierarchy. An entry whose ``needs`` names a key missing here would fail
    at run time and, worse, be absent from the run manifest — so
    ``tests/unit/test_run_manifest_wiring.py`` calls this function to check that
    every registered ontology is resolvable. That guard is only meaningful while
    the test reads the keys from *here* rather than keeping its own copy.
    """
    return {
        "gaf": Path(f"data/raw/goa_annotations/goa_{args.species}.gaf.gz"),
        "go_obo": args.go_ontology,
        "enzyme_dat": args.enzyme_dat,
        "uniprot_dat": args.uniprot_dat,
        "reactome_relations": args.reactome_relations,
        "keywlist": args.keyword_list,
        "subcell": args.subcell,
        "chebi_obo": args.chebi_obo,
        "doid_obo": args.doid_obo,
    }


def start_run_manifest(
    args: argparse.Namespace,
    *,
    ontology_entry: OntologyEntry,
    ontology_label: str,
    ontology_paths: dict,
    interpro_file: Path,
) -> RunManifest:
    """Open the provenance manifest for this run, hashing every input first.

    Which inputs those are comes from what the *selected* registry entry
    declares — ``needs``, plus ``hierarchy_needs`` when True Path propagation is
    on (``src/ontology_registry.py``). Driving it off the registry rather than an
    ``if/elif`` over ``--ontology`` means every registered ontology is covered,
    and a newly registered one is covered without touching this file.

    The manifest is written before the expensive stages, with
    ``"status": "running"``; :meth:`RunManifest.complete` finalizes it. A run
    that fails therefore leaves a visibly unfinished record rather than a stale
    "completed" one from a previous invocation.
    """
    source_urls = input_source_urls(args.species)
    input_records = [
        describe_file(
            interpro_file,
            role="domain_annotations",
            derived_from=source_urls.get("interpro_mappings"),
        )
    ]
    hierarchy_inputs = (
        list(ontology_entry.hierarchy_needs) if args.enable_true_path else []
    )
    # dict.fromkeys de-duplicates while preserving order: an ontology may list
    # the same file as both an annotation and a hierarchy input (subcellular).
    for name in dict.fromkeys(list(ontology_entry.needs) + hierarchy_inputs):
        input_records.append(
            describe_file(
                ontology_paths[name], role=name, source_url=source_urls.get(name)
            )
        )

    return RunManifest(
        args.output_dir / manifest_filename(ontology_label),
        repository=Path.cwd(),
        parameters=vars(args),
        inputs=input_records,
        analysis={
            "ontology": {
                "key": ontology_entry.key,
                "label": ontology_label,
                "ontology_id": ontology_entry.spec.ontology_id,
                "name": ontology_entry.spec.name,
                "term_prefix": ontology_entry.spec.term_prefix,
                "annotation_inputs": list(ontology_entry.needs),
                "supports_true_path": ontology_entry.supports_true_path,
                "true_path_enabled": bool(args.enable_true_path),
                "hierarchy_inputs": hierarchy_inputs,
                "propagation": (
                    "go_dag_with_parental_background_filter"
                    if ontology_entry.external_propagation
                    else "ancestor_closure"
                    if ontology_entry.build_ancestors is not None
                    else None
                ),
            },
            "thresholds": {
                "evidence_filter": args.evidence_filter,
                "fdr_threshold": args.fdr_threshold,
                "fdr_method": "benjamini_hochberg",
                # Single domains and supra-domains are corrected as
                # separate families; each controls FDR at fdr_threshold
                # within itself.
                "fdr_families": ["single", "supra"],
                "fisher_alternative": FISHER_ALTERNATIVE,
                # 0 means no minimum-support filter: an association is kept on
                # FDR significance alone. Recorded explicitly rather than by
                # omission, either way.
                "min_proteins_per_association": args.min_support or None,
                "min_domain_length": MIN_DOMAIN_LENGTH,
                "max_supra_domain_length": MAX_SUPRA_DOMAIN_LENGTH,
                "enable_supra_domains": bool(args.enable_supra_domains),
                "parental_background_min_size": PARENTAL_MIN_BACKGROUND_SIZE,
                "parental_background_alpha": PARENTAL_ALPHA_THRESHOLD,
            },
        },
        command=[sys.executable, *sys.argv],
    )


def main():
    parser = argparse.ArgumentParser(
        description="dcGO Pipeline - Human Protein Analysis",
        # Raw epilog so the ontology table below keeps its line breaks; argument
        # help strings are still wrapped normally.
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="ontologies available to --ontology:\n" + describe_ontologies(),
    )
    parser.add_argument(
        "--domain-key",
        default="interpro",
        choices=list(DOMAIN_KEYS),
        help="Which protein2ipr column defines a domain: 'interpro' (the "
        "integrated entry, default) or 'ssf' (the SUPERFAMILY/SCOP signature, "
        "i.e. the domain universe the published dcGO used). Only the domain "
        "axis changes; the ontology and the statistics are unaffected.",
    )
    parser.add_argument(
        "--evidence-filter",
        default="manual",
        choices=["all", "manual", "experimental"],
        help="GO annotation evidence filter (default: manual)",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.01,
        help="FDR significance threshold (default: 0.01)",
    )
    parser.add_argument(
        "--min-support",
        type=int,
        default=0,
        help="Discard associations supported by fewer than N proteins carrying "
        "both the domain and the term. Applied AFTER the FDR correction, so it "
        "never alters the hypothesis family. Default 0 (no filter): the emergent "
        "domain combinations this method exists to find sit at n = 2-8 proteins, "
        "so a non-zero default would delete them",
    )
    parser.add_argument(
        "--num-cores",
        type=int,
        default=8,
        help="Retained for compatibility with existing scripts and the HPC "
        "batch file. The Fisher stage is compiled Cython and runs in-process, "
        "so this currently has no effect on runtime (default: 8)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Output directory (default: results/)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50000,
        help="Batch size for Fisher tests (default: 50000)",
    )
    parser.add_argument(
        "--species",
        default="human",
        help="Species to analyze: 'human', 'mouse', or specific GOA file name (default: human)",
    )
    parser.add_argument(
        "--ontology",
        default="go",
        choices=ontology_keys(),
        metavar="NAME",
        help="Which ontology to associate domains with (default: go). The "
        "full list, with descriptions, is at the end of this help.",
    )
    parser.add_argument(
        "--xref-db",
        default=None,
        help="UniProt DR database name to harvest when --ontology xref "
        "(e.g. 'KEGG', 'BRENDA', 'GuidetoPHARMACOLOGY')",
    )
    parser.add_argument(
        "--xref-type",
        default=None,
        help="Optional DR third-field filter for --ontology xref "
        "(e.g. 'phenotype' to keep only those typed entries)",
    )
    parser.add_argument(
        "--xref-term-from-type",
        action="store_true",
        help="For --ontology xref, use the DR line's third field as the term "
        "instead of the id (for databases that key the DR line by accession "
        "and carry the vocabulary in that field)",
    )
    parser.add_argument(
        "--enzyme-dat",
        type=Path,
        default=Path("data/raw/enzyme/enzyme.dat"),
        help="Path to Expasy enzyme.dat, used when --ontology ec "
        "(default: data/raw/enzyme/enzyme.dat)",
    )
    parser.add_argument(
        "--uniprot-dat",
        type=Path,
        default=Path("data/raw/uniprot_sprot_dat/uniprot_sprot.dat.gz"),
        help="Path to the UniProt Swiss-Prot flat file, used when --ontology "
        "reactome/keyword/disease/xref "
        "(default: data/raw/uniprot_sprot_dat/uniprot_sprot.dat.gz)",
    )
    parser.add_argument(
        "--enable-true-path",
        action="store_true",
        help="Enable True Path Rule propagation (GO via its OBO DAG, EC via its "
        "numbering, reactome/keyword via their hierarchy files; not available "
        "for disease/xref)",
    )
    parser.add_argument(
        "--go-ontology",
        type=Path,
        default=Path("data/raw/go_ontology/go-basic.obo"),
        help="Path to GO ontology file (default: data/raw/go_ontology/go-basic.obo)",
    )
    parser.add_argument(
        "--reactome-relations",
        type=Path,
        default=Path("data/raw/reactome_relations/ReactomePathwaysRelation.txt"),
        help="Path to Reactome ReactomePathwaysRelation.txt, for --ontology "
        "reactome --enable-true-path",
    )
    parser.add_argument(
        "--keyword-list",
        type=Path,
        default=Path("data/raw/uniprot_keywlist/keywlist.txt"),
        help="Path to UniProt keywlist.txt, for --ontology keyword --enable-true-path",
    )
    parser.add_argument(
        "--subcell",
        type=Path,
        default=Path("data/raw/uniprot_subcell/subcell.txt"),
        help="Path to UniProt subcell.txt (controlled vocabulary + hierarchy), "
        "for --ontology subcellular",
    )
    parser.add_argument(
        "--chebi-obo",
        type=Path,
        default=Path("data/raw/chebi/chebi_lite.obo"),
        help="Path to the ChEBI ontology in OBO format, for --ontology "
        "ligand/cofactor --enable-true-path",
    )
    parser.add_argument(
        "--doid-obo",
        type=Path,
        default=Path("data/raw/disease_ontology/doid.obo"),
        help="Path to the Human Disease Ontology OBO, for --ontology "
        "doid/orphanet_doid (supplies both the OMIM/Orphanet cross-references "
        "used to re-key the annotations and the DAG they propagate up)",
    )
    parser.add_argument(
        "--permute-annotations",
        type=int,
        default=None,
        metavar="SEED",
        help="Calibration control: shuffle which protein carries which term set "
        "(within the analysed universe) before testing. Term and per-protein "
        "annotation marginals are preserved exactly and only the domain↔term "
        "link is broken, so a correctly calibrated run should return ~0 "
        "significant associations. Use it to check that a layer's significant "
        "count reflects signal rather than its hypothesis universe.",
    )
    parser.add_argument(
        "--enable-supra-domains",
        action="store_true",
        default=True,
        help="Include supra-domain (multi-domain) combinations in analysis (default: True)",
    )
    parser.add_argument(
        "--disable-supra-domains",
        dest="enable_supra_domains",
        action="store_false",
        help="Disable supra-domain analysis (single domains only)",
    )

    args = parser.parse_args()

    validate_arguments(args, parser)

    # Validate the arbitrary-cross-reference selection and derive a short label
    # used for logging, the term column, and output filenames. For everything
    # except 'xref' the label is the ontology name (so 'go' stays byte-identical);
    # for 'xref' it is the chosen DR database (e.g. 'kegg').
    if args.ontology == "xref":
        if not args.xref_db:
            parser.error(
                "--ontology xref requires --xref-db (a UniProt DR database name)"
            )
        ontology_label = args.xref_db.lower()
    else:
        ontology_label = args.ontology

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info(f"dcGO PIPELINE - {args.species.upper()} PROTEIN ANALYSIS")
    logger.info("=" * 70)
    logger.info("Configuration:")
    logger.info(f"  Species: {args.species}")
    logger.info(f"  Ontology: {ontology_label.upper()}")
    logger.info(f"  Domain key: {args.domain_key}")
    logger.info(f"  Evidence filter: {args.evidence_filter}")
    logger.info(f"  FDR threshold: {args.fdr_threshold}")
    # Not "CPU cores: N". The Fisher stage is compiled Cython (fisher.pvalue_npy)
    # and runs in-process; --num-cores has no effect on it. Advertising a core
    # count the run does not use is the kind of overstatement the review flagged.
    logger.info(
        f"  CPU cores requested: {args.num_cores} (Fisher stage is single-process)"
    )
    logger.info(
        f"  Supra-domains: {'ENABLED' if args.enable_supra_domains else 'DISABLED'}"
    )
    logger.info(f"  Output directory: {args.output_dir}")

    # Everything the chosen ontology might need, resolved in one place. The
    # registry (src/ontology_registry.py) says which of these it actually uses
    # for its annotations and for its hierarchy.
    ontology_entry = get_ontology(args.ontology)
    ontology_paths = build_ontology_paths(args)

    # True Path Rule propagation needs a term hierarchy: GO's OBO DAG, an
    # implicit one in the term ids (EC, TCDB, MEROPS, CAZy), or a companion
    # hierarchy file (Reactome, keywords, subcellular, ChEBI). Ontologies with
    # no hierarchy (disease, Rhea, xref, …) cannot propagate.
    if args.enable_true_path and not ontology_entry.supports_true_path:
        logger.error(
            f"True Path propagation is not available for --ontology "
            f"{args.ontology} (no term hierarchy). Re-run without "
            "--enable-true-path."
        )
        return 1

    # Fail on missing inputs before the expensive stages rather than degrading
    # silently half-way through a multi-hour run.
    missing = missing_inputs(
        ontology_entry, ontology_paths, for_hierarchy=args.enable_true_path
    )
    if missing:
        logger.error(
            f"--ontology {args.ontology} needs input(s) that are missing: "
            + "; ".join(missing)
        )
        logger.error(
            "Download them: uv run python scripts/download_data.py --list "
            "(then --datasets <name>)"
        )
        return 1

    # File paths - support different species
    interpro_file = Path(f"data/interim/protein2ipr_{args.species}.dat.gz")

    # Check domain annotations exist (shared across ontologies)
    if not interpro_file.exists():
        logger.error(f"{args.species.title()} InterPro file not found: {interpro_file}")
        logger.error(f"Please extract {args.species} data from protein2ipr.dat.gz")
        return 1

    logger.info("Hashing inputs and writing the run provenance manifest...")
    manifest = start_run_manifest(
        args,
        ontology_entry=ontology_entry,
        ontology_label=ontology_label,
        ontology_paths=ontology_paths,
        interpro_file=interpro_file,
    )
    logger.info(f"✓ Run manifest started: {manifest.path}")

    # Build the annotation source for the chosen ontology. Everything downstream
    # only sees a {protein_id: {term}} map, so the engine is ontology-agnostic —
    # see src/annotation_source.py for the seam.
    annotation_source = ontology_entry.build_source(
        ontology_paths,
        {
            "evidence_filter": args.evidence_filter,
            "xref_db": args.xref_db,
            "xref_type": args.xref_type,
            "xref_term_from_type": args.xref_term_from_type,
        },
    )

    # Load data
    logger.info("")
    logger.info("STAGE 1: Loading Data")
    logger.info("─" * 70)

    logger.info(f"Parsing {ontology_label.upper()} annotations...")
    protein_go_map = annotation_source.parse()

    logger.info("Parsing domain annotations...")
    parser_obj = DomainAnnotationParser(
        max_supra_domain_length=MAX_SUPRA_DOMAIN_LENGTH,
        min_domain_length=MIN_DOMAIN_LENGTH,
        domain_key=args.domain_key,
    )
    domain_architectures = parser_obj.parse_protein2ipr_file(interpro_file)

    # Get intersection
    proteins_with_both = set(protein_go_map.keys()) & set(domain_architectures.keys())

    # Restrict the annotation map to that intersection — it defines the protein
    # universe every Fisher table is computed against. See the docstring of
    # restrict_to_universe for why this matters for the UniProt-native sources.
    annotated_proteins = len(protein_go_map)
    protein_go_map = restrict_to_universe(protein_go_map, proteins_with_both)
    if annotated_proteins > len(protein_go_map):
        logger.info(
            f"  Restricted annotations to the domain-annotated universe: "
            f"{annotated_proteins:,} → {len(protein_go_map):,} proteins "
            f"({annotated_proteins - len(protein_go_map):,} dropped, no domain data)"
        )

    # Calibration control. Comparing two ontology layers by their significant
    # counts is only meaningful if neither count is manufactured by its
    # hypothesis universe (test count, term sparsity, marginals). Permuting which
    # protein carries which term set preserves all of that and destroys only the
    # domain↔term relationship, so the significant count under permutation is the
    # layer's own false-positive floor, directly comparable across layers.
    if args.permute_annotations is not None:
        rng = np.random.default_rng(args.permute_annotations)
        proteins = sorted(protein_go_map)
        donors = [proteins[i] for i in rng.permutation(len(proteins))]
        protein_go_map = {
            protein: protein_go_map[donor] for protein, donor in zip(proteins, donors)
        }
        logger.warning(
            "CALIBRATION CONTROL: term sets permuted across "
            f"{len(proteins):,} proteins (seed {args.permute_annotations}). "
            "These results are a null, not predictions."
        )

    # Build protein-domain map (using lists for compatibility with ontology processor)
    # CRITICAL: Include both single domains AND supra-domains as per dcGO methodology
    protein_domain_map = {}
    all_domains = set()
    single_domain_count = 0
    supra_domain_count = 0

    for protein_id in proteins_with_both:
        arch = domain_architectures[protein_id]

        # Always include single domains
        domains = list(arch.single_domains)
        single_domain_count += len(domains)

        # Include supra-domains if enabled (default: True)
        if args.enable_supra_domains:
            domains.extend(arch.supra_domains)
            supra_domain_count += len(arch.supra_domains)

        if domains:
            protein_domain_map[protein_id] = domains
            all_domains.update(domains)

    logger.info(f"  Single domains: {single_domain_count:,} annotations")
    if args.enable_supra_domains:
        logger.info(f"  Supra-domains: {supra_domain_count:,} annotations")
        logger.info(
            f"  Total domain features: {single_domain_count + supra_domain_count:,}"
        )

    # Get all GO terms
    all_go_terms = set()
    for protein_id in proteins_with_both:
        all_go_terms.update(protein_go_map[protein_id])

    domain_list = sorted(all_domains)
    go_list = sorted(all_go_terms)

    logger.info(
        f"✓ Dataset prepared: {len(proteins_with_both):,} proteins, {len(domain_list):,} domains, {len(go_list):,} {ontology_label.upper()} terms"
    )
    logger.info(f"  Total tests: {len(domain_list) * len(go_list):,}")

    # Abort early on an empty design rather than proceeding to zero Fisher tests
    # (which would divide by len(go_list) downstream). Reachable when a
    # species/ontology combination has no overlap between domains and terms.
    if not domain_list or not go_list:
        logger.error(
            f"No domain-{ontology_label.upper()} pairs to test "
            f"({len(domain_list):,} domains, {len(go_list):,} terms). "
            "Check that the annotation and InterPro inputs cover the same proteins."
        )
        return 1

    # Build sparse matrices
    logger.info("")
    logger.info("STAGE 2: Building Sparse Matrices")
    logger.info("─" * 70)
    start_time = time.time()

    protein_domain_matrix, protein_go_matrix, domain_metadata = build_sparse_matrices(
        protein_domain_map, protein_go_map, domain_list, go_list
    )

    matrix_time = time.time() - start_time
    logger.info(f"✓ Sparse matrices built in {matrix_time:.2f}s")

    # Compute contingency tables
    logger.info("")
    logger.info("STAGE 3: Computing Contingency Tables")
    logger.info("─" * 70)
    start_time = time.time()

    tables = compute_contingency_tables_sparse(protein_domain_matrix, protein_go_matrix)

    table_time = time.time() - start_time
    logger.info(
        f"✓ Contingency tables computed in {table_time:.2f}s ({table_time / 60:.1f} min)"
    )

    # Run Fisher's exact tests
    logger.info("")
    logger.info("STAGE 4: Running Fisher's Exact Tests")
    logger.info("─" * 70)
    logger.info(f"Processing {len(tables):,} tests (vectorized Cython, in-process)...")
    start_time = time.time()

    # Progress callback
    def progress_callback(completed, total):
        progress_pct = (completed / total) * 100
        elapsed = time.time() - start_time
        rate = completed / elapsed if elapsed > 0 else 0
        eta = (total - completed) / rate if rate > 0 else 0
        logger.info(
            f"  Progress: {completed:,} / {total:,} ({progress_pct:.1f}%) | {rate:,.0f} tests/s | ETA: {eta / 60:.1f} min"
        )

    odds_ratios, pvalues = fisher_exact_parallel(
        tables,
        alternative=FISHER_ALTERNATIVE,
        n_jobs=args.num_cores,
        batch_size=args.batch_size,
        progress_callback=progress_callback,
    )

    test_time = time.time() - start_time
    logger.info(
        f"✓ Fisher tests completed in {test_time:.2f}s ({test_time / 60:.1f} min)"
    )
    logger.info(f"  Rate: {len(pvalues) / test_time:,.0f} tests/second")

    # Apply FDR correction
    logger.info("")
    logger.info("STAGE 5: FDR Correction")
    logger.info("─" * 70)
    start_time = time.time()

    # Single domains and supra-domains are corrected as separate hypothesis
    # families. A supra-domain is not an exchangeable sibling of its own
    # constituents, and pooling them made the 5.3x larger supra space tighten
    # the threshold for single-domain hypotheses that gain nothing from it.
    # Each family controls FDR at --fdr-threshold within itself.
    is_supra = np.fromiter(
        (
            domain_metadata[domain_id].domain_type is not DomainType.SINGLE
            for domain_id in domain_list
        ),
        dtype=bool,
        count=len(domain_list),
    )
    family = np.where(is_supra, "supra", "single").repeat(len(go_list))

    adjusted_pvalues, thresholds = benjamini_hochberg_by_family(
        pvalues, family, alpha=args.fdr_threshold
    )
    del family, is_supra

    fdr_time = time.time() - start_time
    logger.info(f"✓ FDR correction completed in {fdr_time:.2f}s")
    for label in sorted(thresholds):
        logger.info(f"  Threshold p-value ({label}): {thresholds[label]:.2e}")

    # Count significant associations
    significant = adjusted_pvalues <= args.fdr_threshold

    # Minimum-support filter, applied HERE: after the BH correction, and before
    # anything consumes the significant set.
    #
    # After BH because support is the observed success count that produced the
    # p-value — filtering on it first would shrink the hypothesis family by
    # outcome and leave the q-values anti-conservative, the same mistake fixed
    # in the surprise score in #26. Filtering afterwards only narrows what is
    # reported and changes no q-value.
    #
    # Before the True Path stage because otherwise propagation would run over
    # associations the export then drops, and the propagated file would assert
    # annotations the association file does not support.
    if args.min_support > 0:
        dropped = int((significant & (tables[:, 0, 0] < args.min_support)).sum())
        significant &= tables[:, 0, 0] >= args.min_support
        logger.info(
            f"  Minimum support (n_both >= {args.min_support}, applied after BH): "
            f"{int(significant.sum()):,} kept, {dropped:,} dropped"
        )

    n_significant = int(significant.sum())

    # STAGE 5.5: True Path Rule (Optional)
    propagated_annotations = []
    if args.enable_true_path and n_significant > 0:
        logger.info("")
        logger.info("STAGE 5.5: True Path Rule Propagation")
        logger.info("─" * 70)
        start_time = time.time()

        # Build AssociationResult objects for the significant associations. This
        # is ontology-agnostic — go_list holds GO terms or EC numbers.
        logger.info("Preparing significant associations for propagation...")
        significant_associations = []
        significant_indices = np.where(significant)[0]

        for idx in significant_indices:
            domain_idx = idx // len(go_list)
            go_idx = idx % len(go_list)
            table = tables[idx]
            a, b = int(table[0, 0]), int(table[0, 1])
            c, d = int(table[1, 0]), int(table[1, 1])

            significant_associations.append(
                AssociationResult(
                    domain=domain_list[domain_idx],
                    go_term=go_list[go_idx],
                    p_value=float(pvalues[idx]),
                    q_value=float(adjusted_pvalues[idx]),
                    hyper_score=calculate_hypergeometric_score(a, b, c, d),
                    a=a,
                    b=b,
                    c=c,
                    d=d,
                )
            )

        if ontology_entry.build_ancestors is not None:
            # Every non-GO ontology propagates through the shared engine; only
            # the *ancestors* differ (implicit in EC/TCDB/MEROPS/CAZy ids, or
            # loaded from a hierarchy file for Reactome/keywords/subcellular/
            # ChEBI). See src/ontology_registry.py.
            logger.info(
                f"Propagating {len(significant_associations):,} "
                f"{ontology_label.upper()} associations up the "
                f"{ontology_entry.spec.name} hierarchy..."
            )
            ancestors_fn = ontology_entry.build_ancestors(ontology_paths)
            propagated_annotations = propagate_via_ancestors(
                significant_associations, ancestors_fn
            )
        else:
            # Load GO ontology (OBO DAG) and apply optimal-level filtering.
            logger.info(f"Loading GO ontology from: {args.go_ontology}")
            ontology_processor = OntologyProcessor(args.go_ontology)

            logger.info(
                f"Applying True Path Rule to {len(significant_associations):,} significant associations..."
            )

            # Apply optimal level filtering
            # Note: alpha_threshold is for raw p-values from Fisher tests, not FDR-corrected
            # Using 0.05 as recommended threshold for parent-child comparison tests
            filtered_associations = ontology_processor.apply_optimal_level_filter(
                significant_associations,
                protein_domain_map,
                protein_go_map,
                min_background_size=PARENTAL_MIN_BACKGROUND_SIZE,
                alpha_threshold=PARENTAL_ALPHA_THRESHOLD,
            )

            logger.info(
                f"✓ Optimal level filtering: {len(filtered_associations):,} associations retained"
            )

            # Propagate annotations up the GO hierarchy
            propagated_annotations = ontology_processor.propagate_annotations(
                filtered_associations
            )

        if propagated_annotations:
            direct_count = sum(
                1 for ann in propagated_annotations if ann.annotation_type == "direct"
            )
            propagated_count = len(propagated_annotations) - direct_count

            logger.info(
                f"✓ Generated {len(propagated_annotations):,} total annotations:"
            )
            logger.info(f"  - Direct: {direct_count:,}")
            logger.info(f"  - Propagated: {propagated_count:,}")

        true_path_time = time.time() - start_time
        logger.info(f"✓ True Path Rule completed in {true_path_time:.2f}s")

    # Export results
    logger.info("")
    logger.info("STAGE 6: Exporting Results")
    logger.info("─" * 70)

    # Calculate hypergeometric scores for significant associations
    logger.info("Calculating hypergeometric scores for significant associations...")
    significant_indices = np.where(significant)[0]

    # Output naming is ontology-aware: for GO this reproduces the historical
    # names exactly (domain_go_associations_*.tsv, "go_term" column); EC gets its
    # own files and an "ec_term" column so no consumer of the GO output is affected.
    term_col = f"{ontology_label}_term"
    # Two things qualify the output filename, for the same reason: a run must
    # never silently overwrite results it is meant to be compared against. A
    # non-default domain key distinguishes an SSF-keyed run from the
    # InterPro-keyed one in the same results dir; a permutation seed
    # distinguishes the calibration control from the real run.
    key_prefix = "" if args.domain_key == "interpro" else f"{args.domain_key}_"
    output_label = ontology_label
    if args.permute_annotations is not None:
        output_label = f"{ontology_label}_permuted{args.permute_annotations}"
    assoc_stem = f"domain_{key_prefix}{output_label}_associations"

    # Signature keying is a 1:1 relabelling of InterPro entries over human data,
    # so the entry id comes along for free as a trailing cross-reference column
    # (letting SSF-keyed results be joined to InterPro-keyed ones without a
    # second parse). Nothing is appended for the default keying, so the existing
    # output stays byte-identical.
    xref_header = "" if args.domain_key == "interpro" else "\tinterpro_id"

    def xref_field(domain_id: str) -> str:
        return (
            ""
            if args.domain_key == "interpro"
            else f"\t{parser_obj.interpro_for(domain_id)}"
        )

    # Export significant associations with hypergeometric scores and domain types
    output_file = args.output_dir / f"{assoc_stem}_significant.tsv"
    with open(output_file, "w") as f:
        f.write(
            f"domain\t{term_col}\tp_value\tadj_p_value\todds_ratio\t"
            f"odds_ratio_ci_low\todds_ratio_ci_high\thyper_score\t"
            f"domain_type\tconstituent_domains\tn_observations\ta\tb\tc\td"
            f"{xref_header}\n"
        )
        for idx in significant_indices:
            domain_idx = idx // len(go_list)
            go_idx = idx % len(go_list)
            domain_id = domain_list[domain_idx]

            # Get contingency table values for hypergeometric score
            table = tables[idx]
            a, b = int(table[0, 0]), int(table[0, 1])
            c, d = int(table[1, 0]), int(table[1, 1])
            hyper_score = calculate_hypergeometric_score(a, b, c, d)

            # Get domain metadata
            meta = domain_metadata[domain_id]
            constituents = (
                ",".join(meta.constituent_domains) if meta.constituent_domains else "-"
            )

            or_low, or_high = odds_ratio_interval(a, b, c, d)

            f.write(
                f"{domain_id}\t{go_list[go_idx]}\t"
                f"{pvalues[idx]:.6e}\t{adjusted_pvalues[idx]:.6e}\t{odds_ratios[idx]:.4f}\t"
                f"{or_low:.4f}\t{or_high:.4f}\t{hyper_score:.2f}\t"
                f"{meta.domain_type.value}\t{constituents}\t{meta.observation_count}\t"
                f"{a}\t{b}\t{c}\t{d}"
                f"{xref_field(domain_id)}\n"
            )

    logger.info(f"✓ Exported significant associations to: {output_file}")
    logger.info(f"  {n_significant:,} associations (FDR < {args.fdr_threshold})")

    # Export top associations with hypergeometric scores and domain types
    top_file = args.output_dir / f"{assoc_stem}_top100.tsv"
    top_indices = np.argsort(pvalues)[:100]
    with open(top_file, "w") as f:
        f.write(
            f"rank\tdomain\t{term_col}\tp_value\tadj_p_value\todds_ratio\thyper_score\t"
            f"domain_type\tconstituent_domains\tn_observations{xref_header}\n"
        )
        for rank, idx in enumerate(top_indices, 1):
            domain_idx = idx // len(go_list)
            go_idx = idx % len(go_list)
            domain_id = domain_list[domain_idx]

            # Get contingency table values for hypergeometric score
            table = tables[idx]
            a, b = int(table[0, 0]), int(table[0, 1])
            c, d = int(table[1, 0]), int(table[1, 1])
            hyper_score = calculate_hypergeometric_score(a, b, c, d)

            # Get domain metadata
            meta = domain_metadata[domain_id]
            constituents = (
                ",".join(meta.constituent_domains) if meta.constituent_domains else "-"
            )

            f.write(
                f"{rank}\t{domain_id}\t{go_list[go_idx]}\t"
                f"{pvalues[idx]:.6e}\t{adjusted_pvalues[idx]:.6e}\t{odds_ratios[idx]:.4f}\t{hyper_score:.2f}\t"
                f"{meta.domain_type.value}\t{constituents}\t{meta.observation_count}"
                f"{xref_field(domain_id)}\n"
            )

    logger.info(f"✓ Exported top 100 associations to: {top_file}")

    # Export propagated annotations if True Path Rule was applied
    if propagated_annotations:
        annotations_file = (
            args.output_dir
            / f"domain_{key_prefix}{output_label}_annotations_propagated.tsv"
        )
        with open(annotations_file, "w") as f:
            f.write(
                f"domain\t{term_col}\tq_value\tassociation_score\tannotation_type\tdirect_source_term\n"
            )
            for ann in propagated_annotations:
                f.write(
                    f"{ann.domain}\t{ann.go_term}\t{ann.q_value:.6e}\t{ann.association_score:.2f}\t"
                    f"{ann.annotation_type}\t{ann.direct_source_term}\n"
                )

        logger.info(
            f"✓ Exported {len(propagated_annotations):,} propagated annotations to: {annotations_file}"
        )

    # Performance summary
    total_time = matrix_time + table_time + test_time + fdr_time

    logger.info("")
    logger.info("=" * 70)
    logger.info("PIPELINE COMPLETE!")
    logger.info("=" * 70)
    logger.info("Results Summary:")
    logger.info(f"  Total domain-{ontology_label.upper()} tests: {len(pvalues):,}")
    logger.info(
        f"  Significant associations (FDR < {args.fdr_threshold}): {n_significant:,} ({n_significant / len(pvalues) * 100:.2f}%)"
    )
    if propagated_annotations:
        direct_count = sum(
            1 for ann in propagated_annotations if ann.annotation_type == "direct"
        )
        propagated_count = len(propagated_annotations) - direct_count
        logger.info(
            f"  True Path Rule annotations: {len(propagated_annotations):,} total ({direct_count:,} direct + {propagated_count:,} propagated)"
        )
    logger.info(f"  Total runtime: {total_time:.1f}s ({total_time / 60:.1f} minutes)")
    logger.info("")
    # Finalize the provenance record: output identities plus the counts a reader
    # would otherwise have to trust the log for. Only a run that reaches here is
    # marked "completed" — one that fails or is killed leaves "running".
    output_files = [(output_file, "significant_associations"), (top_file, "top100")]
    if propagated_annotations:
        output_files.append((annotations_file, "propagated_annotations"))
    logger.info("Hashing outputs and finalizing the run manifest...")
    manifest.complete(
        outputs=[describe_file(path, role=role) for path, role in output_files],
        summary={
            "ontology": ontology_label,
            "proteins": len(proteins_with_both),
            "domains": len(domain_list),
            "terms": len(go_list),
            "tests": int(len(pvalues)),
            "significant_associations": n_significant,
            # One cutoff per hypothesis family, not one overall.
            "bh_threshold_pvalue": {k: float(v) for k, v in thresholds.items()},
            "propagated_annotations": len(propagated_annotations),
            "runtime_seconds": round(total_time, 2),
        },
    )

    logger.info("Output files:")
    logger.info(f"  {output_file}")
    logger.info(f"  {top_file}")
    if propagated_annotations:
        logger.info(f"  {annotations_file}")
    logger.info(f"  {manifest.path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
