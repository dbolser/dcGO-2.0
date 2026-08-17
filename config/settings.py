"""
Configuration settings for dcGO Pipeline

This module contains all configuration settings and parameters for the dcGO pipeline,
including data sources, processing parameters, and file paths using Python 3.12 features.
"""

import hashlib
import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Self, Union
from urllib.parse import urlparse

import psutil

from src.release_pins import (
    FLYBASE_FBAL_TO_FBGN_FILENAME,
    FLYBASE_FBGN_UNIPROT_FILENAME,
    FLYBASE_GENOTYPE_PHENOTYPE_FILENAME,
    FLYBASE_RELEASE,
    WORMBASE_PHENOTYPE_FILENAME,
    WORMBASE_RELEASE,
)

logger = logging.getLogger(__name__)


class ConfigurationError(Exception):
    """Exception raised for configuration-related errors."""

    pass


@dataclass(frozen=True)
class DataSource:
    """Configuration for a data source with validation.

    ``checksum`` is an optional ``"<algorithm>:<hex digest>"`` string (a bare
    digest is read as SHA-256). Set it for any source pinned to an immutable
    release URL — ``scripts/download_data.py`` verifies it after downloading, so
    a silently changed or truncated file fails loudly instead of quietly
    altering a run's results. Sources on mutable "current release" URLs leave it
    ``None``, since the file legitimately changes.
    """

    name: str
    url: str
    description: str
    required: bool = True
    #: SHA-256 of the exact bytes we validated against, where the source is a
    #: frozen archive (the published dcGO tables, SCOP 1.75). Left None for
    #: rolling "current_release" URLs, whose content legitimately changes.
    checksum: Optional[str] = None
    #: Size in bytes of those same frozen bytes, so a truncated download or a
    #: silently-substituted file is caught without hashing 90 MB.
    size_bytes: Optional[int] = None
    #: Directory under ``data/raw/`` to save into. Defaults to the source name;
    #: set it when several sources belong together (e.g. the three published
    #: dcGO tables all land in ``data/raw/dcgo_reference/``).
    subdir: Optional[str] = None
    #: What to update when this URL goes stale. Printed by
    #: ``scripts/download_data.py`` on a download failure, for sources whose
    #: URL embeds a release name and therefore *will* break on upstream
    #: release turnover (WormBase, FlyBase).
    update_hint: Optional[str] = None

    def __post_init__(self) -> None:
        """Validate data source configuration."""
        if not self.url:
            raise ConfigurationError(
                f"URL cannot be empty for data source '{self.name}'"
            )

        # Fail on a malformed checksum here, not after a multi-GB download.
        if self.checksum is not None:
            self.checksum_parts()

        # Basic URL validation
        parsed = urlparse(self.url)
        if not parsed.scheme or not parsed.netloc:
            raise ConfigurationError(
                f"Invalid URL format for data source '{self.name}': {self.url}"
            )

        if parsed.scheme not in {"http", "https", "ftp"}:
            raise ConfigurationError(
                f"Unsupported URL scheme for data source '{self.name}': {parsed.scheme}"
            )

    def checksum_parts(self) -> Optional[tuple[str, str]]:
        """``(algorithm, hex digest)`` for :attr:`checksum`, or ``None`` if unset.

        Raises:
            ConfigurationError: the checksum names an algorithm ``hashlib``
                does not provide, or its digest is not hexadecimal.
        """
        if self.checksum is None:
            return None
        algorithm, _, digest = self.checksum.rpartition(":")
        algorithm = (algorithm or "sha256").lower()
        if algorithm not in hashlib.algorithms_available:
            raise ConfigurationError(
                f"Unknown checksum algorithm '{algorithm}' for data source "
                f"'{self.name}'"
            )
        if not digest or any(c not in "0123456789abcdefABCDEF" for c in digest):
            raise ConfigurationError(
                f"Checksum for data source '{self.name}' is not a hex digest: "
                f"{self.checksum!r}"
            )
        return algorithm, digest.lower()


@dataclass(frozen=True)
class ProcessingParameters:
    """Statistical and processing parameters with validation."""

    fdr_threshold: float = 0.01
    #: Not read by the supported ``run_dcgo_human.py`` path, which applies no
    #: minimum-support filter unless ``--min-support`` is given. Kept because
    #: ``DCGO_MIN_PROTEINS`` and the validation below are part of this dataclass'
    #: public surface; see CLAUDE.md for why the default is no filter.
    min_proteins_per_association: int = 3
    max_supra_domain_length: int = 3
    alpha_threshold: float = 0.05
    min_cooccurrence_threshold: int = 5
    max_sequence_length: int = 50000
    min_sequence_length: int = 20
    chunk_size: int = 65536
    timeout: int = 30
    evidence_filter: str = "manual"  # 'all', 'manual' (no IEA), or 'experimental'

    # Supra-domain configuration
    enable_supra_domains: bool = True
    supra_domain_min_count: int = 3
    hierarchical_shrinkage_strength: float = 0.5  # 0-1, higher = more shrinkage
    supra_domain_weighting: str = (
        "empirical_bayes"  # 'none', 'empirical_bayes', 'hierarchical_fdr'
    )

    def __post_init__(self) -> None:
        """Validate processing parameters."""
        if not (0 < self.fdr_threshold < 1):
            raise ConfigurationError(
                f"FDR threshold must be between 0 and 1, got {self.fdr_threshold}"
            )

        if self.min_proteins_per_association <= 0:
            raise ConfigurationError(
                f"Minimum proteins per association must be positive, got {self.min_proteins_per_association}"
            )

        if self.max_supra_domain_length <= 0:
            raise ConfigurationError(
                f"Maximum supra-domain length must be positive, got {self.max_supra_domain_length}"
            )

        if not (0 < self.alpha_threshold < 1):
            raise ConfigurationError(
                f"Alpha threshold must be between 0 and 1, got {self.alpha_threshold}"
            )

        if self.min_cooccurrence_threshold <= 0:
            raise ConfigurationError(
                f"Minimum cooccurrence threshold must be positive, got {self.min_cooccurrence_threshold}"
            )

        if self.max_sequence_length <= self.min_sequence_length:
            raise ConfigurationError(
                f"Maximum sequence length ({self.max_sequence_length}) must be greater than minimum ({self.min_sequence_length})"
            )

        if self.chunk_size <= 0:
            raise ConfigurationError(
                f"Chunk size must be positive, got {self.chunk_size}"
            )

        if self.timeout <= 0:
            raise ConfigurationError(f"Timeout must be positive, got {self.timeout}")

        if self.supra_domain_min_count <= 0:
            raise ConfigurationError(
                f"Minimum supra-domain count must be positive, got {self.supra_domain_min_count}"
            )

        if not (0 <= self.hierarchical_shrinkage_strength <= 1):
            raise ConfigurationError(
                f"Hierarchical shrinkage strength must be between 0 and 1, got {self.hierarchical_shrinkage_strength}"
            )

        if self.supra_domain_weighting not in {
            "none",
            "empirical_bayes",
            "hierarchical_fdr",
        }:
            raise ConfigurationError(
                f"Invalid supra-domain weighting method: {self.supra_domain_weighting}. "
                f"Must be one of: none, empirical_bayes, hierarchical_fdr"
            )


@dataclass(frozen=True)
class ComputeResources:
    """Compute resource configuration with auto-detection and validation."""

    num_cores: int = field(default_factory=lambda: psutil.cpu_count())
    memory_limit_gb: Optional[float] = field(
        default_factory=lambda: psutil.virtual_memory().total / (1024**3)
    )
    java_heap_size: str = "4G"
    temp_dir: Optional[Path] = None

    def __post_init__(self) -> None:
        """Validate compute resource settings."""
        available_cores = psutil.cpu_count()
        if self.num_cores <= 0 or self.num_cores > available_cores:
            raise ConfigurationError(
                f"Number of cores must be between 1 and {available_cores}, got {self.num_cores}"
            )

        if self.memory_limit_gb is not None:
            available_memory = psutil.virtual_memory().total / (1024**3)
            if self.memory_limit_gb <= 0 or self.memory_limit_gb > available_memory:
                raise ConfigurationError(
                    f"Memory limit must be between 0 and {available_memory:.1f}GB, got {self.memory_limit_gb}"
                )

        # Validate Java heap size format
        if not self.java_heap_size.endswith(("M", "G", "m", "g")):
            raise ConfigurationError(
                f"Java heap size must end with 'M' or 'G', got {self.java_heap_size}"
            )

        try:
            int(self.java_heap_size[:-1])
        except ValueError:
            raise ConfigurationError(
                f"Invalid Java heap size format: {self.java_heap_size}"
            )


#: Species names that mean "the multi-species background", not an organism.
#: These resolve to EBI's single cross-organism GOA release rather than to a
#: per-species directory — see ``Config.goa_url_for``.
ALL_SPECIES_ALIASES = frozenset({"allspecies", "all_species", "uniprot_all"})

#: Species names of the form ``<base>_t0_2021`` / ``<base>_t1_2023``: local
#: snapshots pinned to a numbered GOA archive release for the temporal
#: benchmark. The release number is not recoverable from the name, so no
#: upstream URL can be composed for them — see ``Config.goa_url_for``.
DERIVED_SNAPSHOT_SPECIES = re.compile(r"_t[01]_\d{4}$")


@dataclass
class Config:
    """
    Main configuration class for dcGO Pipeline using Python 3.12 dataclass features.

    This class provides comprehensive configuration management with:
    - Type-safe configuration parameters
    - Automatic validation
    - Environment variable override support
    - Path management with automatic directory creation
    - Resource auto-detection
    """

    # Core paths (calculated from current file location)
    base_dir: Path = field(
        default_factory=lambda: Path(__file__).parent.parent.resolve()
    )

    # Base URL for dated (archived) human GOA snapshots, one numbered release per
    # file (e.g. goa_human.gaf.205.gz). Used by scripts/download_data.py
    # --goa-archive for the temporal benchmark (VALIDATION_PLAN §2/§6).
    goa_archive_base_url: str = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN"

    # Base URL for current per-species GOA releases. EBI lays these out as
    # <base>/<SPECIES_UPPER>/goa_<species>.gaf.gz (HUMAN/goa_human.gaf.gz,
    # MOUSE/goa_mouse.gaf.gz, ZEBRAFISH/goa_zebrafish.gaf.gz, …). Used by
    # scripts/download_data.py --species to fetch non-human annotations.
    goa_base_url: str = "https://ftp.ebi.ac.uk/pub/databases/GO/goa"

    # Data sources configuration
    data_sources: Dict[str, DataSource] = field(
        default_factory=lambda: {
            "uniprot_sprot": DataSource(
                name="uniprot_sprot",
                url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
                description="UniProt Swiss-Prot protein sequences (only needed for local HMM scanning)",
                required=False,
            ),
            "uniprot_trembl": DataSource(
                name="uniprot_trembl",
                url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
                description="UniProt TrEMBL protein sequences",
                required=False,
            ),
            # UniProt Swiss-Prot flat file (DR cross-references + KW keywords).
            # UniProt is the protein universe, so this is the source of
            # UniProt-native term annotations that need no identifier mapping:
            # Reactome, KEGG, keywords, disease DBs, etc. Consumed by
            # src/uniprot_annotation_source.py (run_dcgo_human.py --ontology
            # reactome|keyword). Large (~1 GB compressed).
            "uniprot_sprot_dat": DataSource(
                name="uniprot_sprot_dat",
                url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
                description="UniProt Swiss-Prot flat file (DR cross-refs + keywords) for UniProt-native ontologies",
                required=False,
            ),
            # Reactome pathway hierarchy (parent<TAB>child), for True Path
            # propagation of --ontology reactome. Consumed by
            # src/uniprot_annotation_source.parse_reactome_relations.
            "reactome_relations": DataSource(
                name="reactome_relations",
                url="https://reactome.org/download/current/ReactomePathwaysRelation.txt",
                description="Reactome pathway parent/child relations (True Path for --ontology reactome)",
                required=False,
            ),
            # UniProt keyword list (defines the keyword hierarchy on HI lines),
            # for True Path propagation of --ontology keyword. Consumed by
            # src/uniprot_annotation_source.parse_keyword_hierarchy.
            "uniprot_keywlist": DataSource(
                name="uniprot_keywlist",
                url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/keywlist.txt",
                description="UniProt keyword list + hierarchy (True Path for --ontology keyword)",
                required=False,
            ),
            # UniProt subcellular-location controlled vocabulary: maps the
            # location strings in CC SUBCELLULAR LOCATION comments to stable
            # SL- accessions and defines their hierarchy (HI/HP lines).
            # Consumed by src/uniprot_annotation_source.parse_subcell_vocabulary
            # (run_dcgo_human.py --ontology subcellular).
            "uniprot_subcell": DataSource(
                name="uniprot_subcell",
                url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/subcell.txt",
                description="UniProt subcellular location vocabulary + hierarchy (--ontology subcellular)",
                required=False,
            ),
            # ChEBI chemical ontology, for True Path propagation of the ligand
            # and cofactor layers (--ontology ligand|cofactor). The "lite"
            # flavour carries the is_a graph without the chemical data blocks.
            "chebi": DataSource(
                name="chebi",
                url="https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo",
                description="ChEBI ontology (True Path for --ontology ligand|cofactor)",
                required=False,
            ),
            # Human Disease Ontology. Supplies both the DOID is_a DAG (True Path
            # for --ontology doid|orphanet_doid) and the OMIM/Orphanet
            # cross-references that re-key UniProt's disease layer onto it.
            # Consumed by src/disease_ontology.py.
            #
            # Pinned to an immutable OBO Foundry *release* PURL rather than the
            # mutable https://purl.obolibrary.org/obo/doid.obo, so a run is
            # reproducible; the checksum is verified on download. Bump both
            # together when refreshing (the release date is the OBO header's
            # data-version).
            "disease_ontology": DataSource(
                name="disease_ontology",
                url="https://purl.obolibrary.org/obo/doid/releases/2026-07-31/doid.obo",
                description="Human Disease Ontology (DOID DAG + OMIM/Orphanet xrefs) for --ontology doid",
                required=False,
                checksum="sha256:5b9803aa17eeabf4c70f144c64216294d01e66335da3e560576d4eb2dc9ff490",
            ),
            # Mondo Disease Ontology. Like doid.obo it supplies both the term
            # DAG and the OMIM/Orphanet cross-references (prefixes OMIM: and
            # Orphanet:, vs DO's MIM:/ORDO:) that re-key UniProt's disease
            # layer, for --ontology mondo|orphanet_mondo. Pinned to an
            # immutable release PURL; bump URL and checksum together.
            "mondo_ontology": DataSource(
                name="mondo_ontology",
                url="https://purl.obolibrary.org/obo/mondo/releases/2026-08-04/mondo.obo",
                description="Mondo Disease Ontology (DAG + OMIM/Orphanet xrefs) for --ontology mondo",
                required=False,
                subdir="mondo",
                checksum="sha256:fce4bccd97c4eb66161e88272ca0bd6ecd28c003afa878bbec8eae0ceb78fba8",
            ),
            # Human Phenotype Ontology annotation and ontology files, for
            # run_dcgo_human.py --ontology hpo (src/hpo_annotation_source.py).
            # Both land in data/raw/hpo/. genes_to_phenotype.txt is NCBI-GeneID
            # keyed; the GeneID → accession re-key uses the Swiss-Prot flat
            # file (uniprot_sprot_dat above), so no idmapping download.
            "hpo_annotations": DataSource(
                name="hpo_annotations",
                url="https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt",
                description="HPO gene → phenotype annotations (--ontology hpo)",
                required=False,
                subdir="hpo",
            ),
            "hpo_ontology": DataSource(
                name="hpo_ontology",
                url="https://purl.obolibrary.org/obo/hp.obo",
                description="Human Phenotype Ontology DAG (True Path for --ontology hpo)",
                required=False,
                subdir="hpo",
            ),
            # GWAS Catalog association file, EFO-annotated, for --ontology efo
            # (src/gwas_annotation_source.py). The per-file .tsv names on the
            # FTP are gone; associations now ship as one zip under
            # releases/latest/ (a moving target — the run manifest's SHA-256
            # identifies the actual snapshot; 2026-08-02 at acquisition).
            "gwas_catalog": DataSource(
                name="gwas_catalog",
                url="https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
                "gwas-catalog-associations_ontology-annotated-full.zip",
                description="GWAS Catalog SNP → EFO trait associations (--ontology efo)",
                required=False,
                subdir="gwas_catalog",
            ),
            # EFO ships versioned releases on GitHub; pinned like the other
            # ontology releases (bump URL and checksum together).
            "efo_ontology": DataSource(
                name="efo_ontology",
                url="https://github.com/EBISPOT/efo/releases/download/v3.93.0/efo.obo",
                description="Experimental Factor Ontology DAG (True Path for --ontology efo)",
                required=False,
                subdir="efo",
                checksum="sha256:66e87fc65a6254c6d69281ed3d286784ee5f8265b1e57691efddd29b20570c46",
            ),
            # HPA single-cell expression matrix, for --ontology celltype
            # (src/hpa_annotation_source.py). Unversioned "current release"
            # URL (HPA v25 at acquisition, 2025-12-12); the run manifest's
            # SHA-256 identifies the snapshot.
            "hpa_single_cell": DataSource(
                name="hpa_single_cell",
                url="https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip",
                description="HPA single-cell gene × cell-type expression (--ontology celltype)",
                required=False,
                subdir="hpa",
            ),
            # SynGO bulk release: one zip carrying both the HGNC-keyed
            # annotations and the term hierarchy, for --ontology syngo
            # (src/syngo_annotation_source.py). Release 1.3 (2025-03) is the
            # current bulk download; the portal itself is a JS app.
            "syngo": DataSource(
                name="syngo",
                url="https://www.syngoportal.org/data/syngo1.3_complete_data.zip",
                description="SynGO synaptic annotations + ontology (--ontology syngo)",
                required=False,
            ),
            # ---- Model-organism phenotype layers ---------------------------
            # Each pairs an annotation file (gene → phenotype term) with the
            # database's own gene → UniProt id-mapping file, and an OBO for
            # True Path / relative inference. They run against the organism's
            # protein universe: --species mouse/worm/zebrafish/fly.
            "mgi_genepheno": DataSource(
                name="mgi_genepheno",
                url="https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt",
                description="MGI mouse genotype → MP phenotype report (--ontology mp)",
                required=False,
                subdir="mgi_reports",
            ),
            "mgi_marker_swissprot": DataSource(
                name="mgi_marker_swissprot",
                url="https://www.informatics.jax.org/downloads/reports/MRK_SwissProt_TrEMBL.rpt",
                description="MGI marker → UniProt accession report (--ontology mp)",
                required=False,
                subdir="mgi_reports",
            ),
            "mp_ontology": DataSource(
                name="mp_ontology",
                url="https://purl.obolibrary.org/obo/mp.obo",
                description="Mammalian Phenotype Ontology DAG (True Path for --ontology mp)",
                required=False,
            ),
            # The filename carries the WormBase release, and the URL WILL 404
            # when WormBase moves to the next release: only the
            # current-production-release/ directory is servable (the
            # per-release archive paths 403, and no unversioned filename is
            # published — verified against the directory listing). The release
            # string is pinned once in src/release_pins.py, shared with the
            # run_dcgo_human.py flag default, so a bump is one edit there.
            # The run manifest labels inputs with this URL; after a bump the
            # file's SHA-256 remains its identity.
            "wormbase_phenotype": DataSource(
                name="wormbase_phenotype",
                url="https://downloads.wormbase.org/releases/"
                "current-production-release/ONTOLOGY/"
                f"{WORMBASE_PHENOTYPE_FILENAME}",
                description="WormBase gene → phenotype GAF (--ontology wbphenotype)",
                required=False,
                subdir="wormbase",
                update_hint=(
                    f"WormBase has likely moved past {WORMBASE_RELEASE}: bump "
                    "WORMBASE_RELEASE in src/release_pins.py (updates this URL "
                    "and the --wb-phenotype default together)"
                ),
            ),
            "wbphenotype_ontology": DataSource(
                name="wbphenotype_ontology",
                url="https://purl.obolibrary.org/obo/wbphenotype.obo",
                description="WormBase Phenotype Ontology DAG (True Path for --ontology wbphenotype)",
                required=False,
                subdir="wormbase_ontology",
            ),
            # UniProt's per-organism idmapping covers TrEMBL too, which the
            # Swiss-Prot flat file route used by the human gene-keyed layers
            # does not — and model-organism proteomes are mostly TrEMBL.
            "worm_idmapping": DataSource(
                name="worm_idmapping",
                url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/CAEEL_6239_idmapping.dat.gz",
                description="UniProt idmapping for C. elegans (WBGene → accession)",
                required=False,
                subdir="uniprot_idmapping",
            ),
            # The other model organisms' idmapping files. The mp/zfa/fbcv/fbbt
            # annotation chains map through their database's own tables, but
            # scripts/extract_species_interpro.py builds each species' domain
            # universe from these, so the documented per-species chain needs
            # them downloadable too.
            "mouse_idmapping": DataSource(
                name="mouse_idmapping",
                url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz",
                description="UniProt idmapping for mouse (accession universe for --species mouse)",
                required=False,
                subdir="uniprot_idmapping",
            ),
            "zebrafish_idmapping": DataSource(
                name="zebrafish_idmapping",
                url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/DANRE_7955_idmapping.dat.gz",
                description="UniProt idmapping for zebrafish (accession universe for --species zebrafish)",
                required=False,
                subdir="uniprot_idmapping",
            ),
            "fly_idmapping": DataSource(
                name="fly_idmapping",
                url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/DROME_7227_idmapping.dat.gz",
                description="UniProt idmapping for D. melanogaster (accession universe for --species fly)",
                required=False,
                subdir="uniprot_idmapping",
            ),
            "zfin_phenotype": DataSource(
                name="zfin_phenotype",
                url="https://zfin.org/downloads/phenoGeneCleanData_fish.txt",
                description="ZFIN clean gene → EQ phenotype data (--ontology zfa)",
                required=False,
                subdir="zfin",
            ),
            "zfin_uniprot": DataSource(
                name="zfin_uniprot",
                url="https://zfin.org/downloads/uniprot.txt",
                description="ZFIN ZDB-GENE → UniProt accession mapping (--ontology zfa)",
                required=False,
                subdir="zfin",
            ),
            "zfa_ontology": DataSource(
                name="zfa_ontology",
                url="https://purl.obolibrary.org/obo/zfa.obo",
                description="Zebrafish Anatomy Ontology DAG (True Path for --ontology zfa)",
                required=False,
                subdir="zfin_ontology",
            ),
            # FlyBase paths and filenames carry the release; pinned once in
            # src/release_pins.py (shared with the run_dcgo_human.py flag
            # defaults, so a release bump is one edit there). Unlike WormBase,
            # FlyBase keeps old releases servable, so these URLs stay valid
            # after turnover. The working host is s3ftp.flybase.org
            # (ftp.flybase.net TLS-fails). The run manifest labels inputs
            # with these URLs; the file's SHA-256 is its identity.
            "flybase_genotype_phenotype": DataSource(
                name="flybase_genotype_phenotype",
                url=f"https://s3ftp.flybase.org/releases/{FLYBASE_RELEASE}/"
                "precomputed_files/alleles/"
                f"{FLYBASE_GENOTYPE_PHENOTYPE_FILENAME}",
                description="FlyBase genotype → FBcv/FBbt phenotype table (--ontology fbcv/fbbt)",
                required=False,
                subdir="flybase",
                update_hint=(
                    "bump FLYBASE_RELEASE in src/release_pins.py (updates the "
                    "three FlyBase URLs and their flag defaults together)"
                ),
            ),
            "flybase_fbal_to_fbgn": DataSource(
                name="flybase_fbal_to_fbgn",
                url=f"https://s3ftp.flybase.org/releases/{FLYBASE_RELEASE}/"
                f"precomputed_files/alleles/{FLYBASE_FBAL_TO_FBGN_FILENAME}",
                description="FlyBase allele → gene mapping (--ontology fbcv/fbbt)",
                required=False,
                subdir="flybase",
                update_hint="bump FLYBASE_RELEASE in src/release_pins.py",
            ),
            "flybase_fbgn_uniprot": DataSource(
                name="flybase_fbgn_uniprot",
                url=f"https://s3ftp.flybase.org/releases/{FLYBASE_RELEASE}/"
                f"precomputed_files/genes/{FLYBASE_FBGN_UNIPROT_FILENAME}",
                description="FlyBase FBgn → UniProt accession mapping (--ontology fbcv/fbbt)",
                required=False,
                subdir="flybase",
                update_hint="bump FLYBASE_RELEASE in src/release_pins.py",
            ),
            "fbbt_ontology": DataSource(
                name="fbbt_ontology",
                url="https://purl.obolibrary.org/obo/fbbt.obo",
                description="Drosophila Anatomy Ontology DAG (True Path for --ontology fbbt)",
                required=False,
                subdir="flybase_ontology",
            ),
            "fbcv_ontology": DataSource(
                name="fbcv_ontology",
                url="https://purl.obolibrary.org/obo/fbcv.obo",
                description="FlyBase Controlled Vocabulary DAG (True Path for --ontology fbcv)",
                required=False,
                subdir="flybase_ontology",
            ),
            "goa_annotations": DataSource(
                name="goa_annotations",
                url="https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz",
                description="Gene Ontology Annotation (GOA) database for human",
                required=True,
            ),
            "go_ontology": DataSource(
                name="go_ontology",
                url="https://current.geneontology.org/ontology/go-basic.obo",
                description="Gene Ontology basic ontology file",
                required=True,
            ),
            # Pre-computed domain annotations (preferred approach)
            "interpro_mappings": DataSource(
                name="interpro_mappings",
                url="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz",
                description="Pre-computed InterPro domain mappings for UniProt proteins",
                required=True,
            ),
            "pfam_regions": DataSource(
                name="pfam_regions",
                url="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.regions.tsv.gz",
                description="Pfam domain region annotations",
                required=False,
            ),
            "interpro_definitions": DataSource(
                name="interpro_definitions",
                url="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro.xml.gz",
                description="InterPro entry definitions and metadata",
                required=False,
            ),
            # Expasy ENZYME database for the EC ontology path. Maps each EC
            # number to its UniProt accessions (DR lines) — same id space as
            # protein2ipr, so no identifier mapping is needed. Consumed by
            # src/ec_annotation_source.py (run_dcgo_human.py --ontology ec).
            "enzyme": DataSource(
                name="enzyme",
                url="https://ftp.expasy.org/databases/enzyme/enzyme.dat",
                description="Expasy ENZYME (EC number → UniProt accession) for the EC ontology path",
                required=False,
            ),
            # Curated InterPro->GO mapping used as a validation reference (§1)
            "interpro2go": DataSource(
                name="interpro2go",
                url="https://current.geneontology.org/ontology/external2go/interpro2go",
                description="Manually curated InterPro2GO mappings (validation reference)",
                required=False,
            ),
            # ---- Published dcGO (Fang & Gough 2013) reference tables --------
            # The comparator for VALIDATION_PLAN §3. Still served by SUPERFAMILY
            # at supfam.org; the /SUPERFAMILY/ path prefix is required (the
            # shorter https://supfam.org/cgi-bin/dcdownload.cgi returns 500), and
            # the human-facing download index (cgi-bin/dcdownload.cgi) takes
            # ~36 s to render — these Domain2GO/ files themselves are fast.
            # Domains are keyed by bare SCOP sunid, which is what our
            # --domain-key ssf run produces (SSFnnnnn → sunid nnnnn).
            # Consumed by validation/compare_original_dcgo.py.
            "dcgo_domain2go_sql": DataSource(
                name="dcgo_domain2go_sql",
                url="https://supfam.org/SUPERFAMILY/Domain2GO/Domain2GO.sql.gz",
                description="Published dcGO domain→GO MySQL dump (GO_mapping carries per-pair FDR + h-score)",
                required=False,
                checksum="9681c5a68cf35d48ca2a2de1921560aacee861694bb87ae4bbb881fe24db714b",
                size_bytes=15810360,
                subdir="dcgo_reference",
            ),
            "dcgo_sp2go": DataSource(
                name="dcgo_sp2go",
                url="https://supfam.org/SUPERFAMILY/Domain2GO/SP2GO.txt",
                description="Published dcGO supra-domain→GO associations (comma-joined sunids)",
                required=False,
                checksum="cdb2d6f073ec53c006091fd32b2f0c7eaf1e16ea1f20df3975603abea5839862",
                size_bytes=92007833,
                subdir="dcgo_reference",
            ),
            "dcgo_domain2go_flat": DataSource(
                name="dcgo_domain2go_flat",
                url="https://supfam.org/SUPERFAMILY/Domain2GO/Domain2GO_supported_only_by_all.txt",
                description="Published dcGO domain→GO flat file (high-coverage 'all proteins' set, with IC and direct/inherited flag)",
                required=False,
                checksum="534141679c4d85705bb94ed8d91d9829cb8f8bd51c4fe4d9b6c3643936d69a5c",
                size_bytes=37270998,
                subdir="dcgo_reference",
            ),
            # ---- SCOP 1.75 -------------------------------------------------
            # InterPro's SUPERFAMILY member database is pinned to the SCOP 1.75
            # HMM library — the same release the 2013 dcGO used — so 1.75 (not
            # SCOPe 2.08) is the right release for resolving our sunids. Note the
            # legacy naming: dir.<x>.scop.1.75.txt (dir.<x>.scop.txt_1.75 404s).
            "scop_des": DataSource(
                name="scop_des",
                url="https://scop.berkeley.edu/downloads/parse/dir.des.scop.1.75.txt",
                description="SCOP 1.75 node descriptions (sunid → type/sccs/name), for SSF identifier resolution",
                required=False,
                checksum="1f90e45bd527a433c938acb619ae83f61837d1435e5b71fc868fa43bc2d6c18c",
                size_bytes=6029837,
                subdir="scop",
            ),
            "scop_hie": DataSource(
                name="scop_hie",
                url="https://scop.berkeley.edu/downloads/parse/dir.hie.scop.1.75.txt",
                description="SCOP 1.75 hierarchy (sunid → parent → children), class > fold > superfamily > family",
                required=False,
                checksum="59f763ae61b27757eb8e966270e9a2ddaa8c9ed7d6a9b3a94e0186b49036804e",
                size_bytes=2961802,
                subdir="scop",
            ),
            # Optional: Local computation tools
            "pfam_hmms": DataSource(
                name="pfam_hmms",
                url="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
                description="Pfam HMM profiles for local domain detection",
                required=False,
            ),
            "interpro_scan": DataSource(
                name="interpro_scan",
                url="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/interproscan-5.67-99.0-64-bit.tar.gz",
                description="InterProScan software for local domain analysis",
                required=False,
            ),
        }
    )

    # Processing parameters
    processing: ProcessingParameters = field(default_factory=ProcessingParameters)

    # Compute resources
    compute: ComputeResources = field(default_factory=ComputeResources)

    # Environment overrides
    use_env_overrides: bool = True

    def __post_init__(self) -> None:
        """Initialize configuration with validation and environment overrides."""
        # Apply environment variable overrides if enabled
        if self.use_env_overrides:
            self._apply_env_overrides()

        # Create directory structure
        self._create_directories()

        # Log configuration summary
        self._log_configuration()

    @property
    def data_dir(self) -> Path:
        """Directory for storing raw and processed data files."""
        return self.base_dir / "data"

    @property
    def results_dir(self) -> Path:
        """Directory for storing analysis results."""
        return self.base_dir / "results"

    @property
    def logs_dir(self) -> Path:
        """Directory for storing log files."""
        return self.base_dir / "logs"

    @property
    def cache_dir(self) -> Path:
        """Directory for storing cached intermediate results."""
        return self.data_dir / "cache"

    @property
    def database_path(self) -> Path:
        """Path to the main SQLite database file."""
        return self.results_dir / "dcgo_database.db"

    @property
    def temp_dir(self) -> Path:
        """Temporary directory for processing operations."""
        if self.compute.temp_dir:
            return self.compute.temp_dir
        return self.base_dir / "temp"

    def _apply_env_overrides(self) -> None:
        """Apply environment variable overrides to configuration."""
        # Processing parameter overrides
        if fdr_env := os.getenv("DCGO_FDR_THRESHOLD"):
            try:
                fdr_value = float(fdr_env)
                object.__setattr__(self.processing, "fdr_threshold", fdr_value)
            except ValueError:
                logger.warning(f"Invalid FDR threshold in environment: {fdr_env}")

        if cores_env := os.getenv("DCGO_NUM_CORES"):
            try:
                cores_value = int(cores_env)
                object.__setattr__(self.compute, "num_cores", cores_value)
            except ValueError:
                logger.warning(f"Invalid number of cores in environment: {cores_env}")

        if min_proteins_env := os.getenv("DCGO_MIN_PROTEINS"):
            try:
                min_proteins_value = int(min_proteins_env)
                object.__setattr__(
                    self.processing, "min_proteins_per_association", min_proteins_value
                )
            except ValueError:
                logger.warning(
                    f"Invalid minimum proteins in environment: {min_proteins_env}"
                )

        if java_heap_env := os.getenv("DCGO_JAVA_HEAP"):
            object.__setattr__(self.compute, "java_heap_size", java_heap_env)

    def _create_directories(self) -> None:
        """Create necessary directories if they don't exist."""
        directories = [
            self.data_dir / "raw",
            self.data_dir / "processed",
            self.data_dir / "interim",
            self.cache_dir,
            self.results_dir,
            self.logs_dir,
            self.temp_dir,
        ]

        for directory in directories:
            try:
                directory.mkdir(parents=True, exist_ok=True)
                logger.debug(f"Ensured directory exists: {directory}")
            except PermissionError as e:
                raise ConfigurationError(f"Cannot create directory {directory}: {e}")

    def _log_configuration(self) -> None:
        """Log configuration summary for debugging."""
        logger.info("dcGO Pipeline Configuration Summary:")
        logger.info(f"  Base directory: {self.base_dir}")
        logger.info(f"  Data directory: {self.data_dir}")
        logger.info(f"  Results directory: {self.results_dir}")
        logger.info(f"  CPU cores: {self.compute.num_cores}")
        logger.info(f"  Memory limit: {self.compute.memory_limit_gb:.1f}GB")
        logger.info(f"  FDR threshold: {self.processing.fdr_threshold}")
        logger.info(
            f"  Required data sources: {len([ds for ds in self.data_sources.values() if ds.required])}"
        )

    @classmethod
    def from_dict(cls, config_dict: Dict) -> Self:
        """Create configuration from dictionary."""
        return cls(**config_dict)

    @classmethod
    def from_env(cls) -> Self:
        """Create configuration with environment variable overrides enabled."""
        return cls(use_env_overrides=True)

    def goa_url_for(self, species: str) -> str:
        """Build the current-release GOA download URL for a species.

        EBI publishes per-species GOA under ``<base>/<SPECIES_UPPER>/`` with a
        ``goa_<species>.gaf.gz`` filename, e.g. ``MOUSE/goa_mouse.gaf.gz``. For
        ``human`` this reproduces the URL pinned in ``data_sources``.

        Two kinds of name are not organisms and must not be run through that
        pattern, because it would compose a plausible-looking URL that 404s and
        every run manifest would then record it as the input's origin:

        * the multi-species background (``allspecies`` and friends), whose
          upstream is the single cross-organism release
          ``UNIPROT/goa_uniprot_all.gaf.gz``;
        * temporal snapshots such as ``human_t0_2021`` or
          ``allspecies_t0_2021``, which are local files pinned to a *numbered*
          archive release (``goa_human.gaf.205.gz``). The release number is not
          recoverable from the species name, so there is no URL to give and
          this raises rather than invent one.

        Raises:
            ConfigurationError: if the species is empty, or is a temporal
                snapshot with no composable upstream URL.
        """
        species = species.strip().lower()
        if not species:
            raise ConfigurationError("Species must be a non-empty string")
        if DERIVED_SNAPSHOT_SPECIES.search(species):
            raise ConfigurationError(
                f"{species!r} is a temporal snapshot pinned to a numbered GOA "
                f"archive release, not a species directory; its upstream URL "
                f"cannot be composed from the name. Fetch it with "
                f"--goa-archive and record the release explicitly."
            )
        if species in ALL_SPECIES_ALIASES:
            return f"{self.goa_base_url}/UNIPROT/goa_uniprot_all.gaf.gz"
        return f"{self.goa_base_url}/{species.upper()}/goa_{species}.gaf.gz"

    def get_data_source_url(self, source_name: str) -> str:
        """Get URL for a specific data source."""
        if source_name not in self.data_sources:
            raise ConfigurationError(f"Unknown data source: {source_name}")
        return self.data_sources[source_name].url

    def get_required_data_sources(self) -> List[str]:
        """Get list of required data source names."""
        return [name for name, source in self.data_sources.items() if source.required]

    def validate_paths(self) -> bool:
        """Validate all configured paths are accessible."""
        paths_to_check = [self.base_dir, self.data_dir, self.results_dir, self.logs_dir]

        for path in paths_to_check:
            if not path.exists():
                logger.error(f"Required path does not exist: {path}")
                return False

            if not os.access(path, os.R_OK | os.W_OK):
                logger.error(f"Insufficient permissions for path: {path}")
                return False

        return True

    def get_system_info(self) -> Dict[str, Union[str, int, float]]:
        """Get system information for diagnostics."""
        return {
            "python_version": os.sys.version,
            "cpu_cores_available": psutil.cpu_count(),
            "cpu_cores_configured": self.compute.num_cores,
            "memory_total_gb": psutil.virtual_memory().total / (1024**3),
            "memory_configured_gb": self.compute.memory_limit_gb,
            "disk_free_gb": psutil.disk_usage(self.base_dir).free / (1024**3),
            "base_directory": str(self.base_dir),
        }

    # Legacy compatibility properties for existing code
    @property
    def DATASOURCES(self) -> Dict[str, str]:
        """Legacy compatibility: return data sources as URL dictionary."""
        return {name: source.url for name, source in self.data_sources.items()}

    @property
    def FDR_THRESHOLD(self) -> float:
        """Legacy compatibility: FDR threshold."""
        return self.processing.fdr_threshold

    @property
    def MIN_PROTEINS_PER_ASSOCIATION(self) -> int:
        """Legacy compatibility: minimum proteins per association."""
        return self.processing.min_proteins_per_association

    @property
    def MAX_SUPRA_DOMAIN_LENGTH(self) -> int:
        """Legacy compatibility: maximum supra-domain length."""
        return self.processing.max_supra_domain_length

    @property
    def NUM_CORES(self) -> int:
        """Legacy compatibility: number of CPU cores."""
        return self.compute.num_cores

    @property
    def BASE_DIR(self) -> Path:
        """Legacy compatibility: base directory."""
        return self.base_dir

    @property
    def DATA_DIR(self) -> Path:
        """Legacy compatibility: data directory."""
        return self.data_dir

    @property
    def RESULTS_DIR(self) -> Path:
        """Legacy compatibility: results directory."""
        return self.results_dir

    @property
    def LOGS_DIR(self) -> Path:
        """Legacy compatibility: logs directory."""
        return self.logs_dir

    @property
    def DATABASE_PATH(self) -> Path:
        """Legacy compatibility: database path."""
        return self.database_path


# Global configuration instance for backward compatibility
config = Config()
