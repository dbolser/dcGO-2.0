"""Generic id-remapping for annotation maps, with audited coverage.

Several layers translate one id space into another *before* the statistics see
the annotations: the DOID layer re-keys OMIM/Orphanet disease ids onto Disease
Ontology terms (the term axis), the HPO and SynGO layers re-key NCBI/HGNC gene
ids onto UniProt accessions (the protein axis). The mechanics and the audit
are identical — map each value through a translation table, expand one-to-many
targets, drop-and-count unmapped ids — so they live here, neutrally named,
rather than in any one ontology's module.

Mapping coverage is the first-class number for every such layer: a re-key is
only worth having if it reaches most of the annotations, and reporting "more
significant associations" without it would be meaningless, because dropping
unmapped ids shrinks the hypothesis universe on its own.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Protocol, Set, Tuple

from loguru import logger


class SupportsTargets(Protocol):
    """A source-id → target-ids translation table, as :func:`remap_values` sees it.

    :class:`src.disease_ontology.XrefMapping` is the DOID implementation; the
    gene→accession maps in :mod:`src.gene_mapping` are another.
    """

    def targets(self, source_id: str) -> Set[str]: ...


@dataclass
class RemapCoverage:
    """What re-keying one axis of a ``{key: {value}}`` annotation map did.

    The *values* are the ids being remapped (OMIM ids for the DOID layer, gene
    ids for the gene-keyed layers); the *keys* are the untouched other axis
    (proteins for DOID, ontology terms for the gene-keyed layers, which remap
    the inverted map).

    Attributes:
        n_source_values: distinct source ids in the input map.
        n_mapped_values: of those, how many had at least one target.
        n_source_annotations: key→value pairs in the input map.
        n_mapped_annotations: input pairs whose value mapped (the coverage that
            actually matters — an id used by many keys counts many times).
        n_result_values: distinct target ids in the output map.
        n_result_annotations: key→target pairs in the output (can exceed
            ``n_mapped_annotations`` through one-to-many expansion, or fall
            below it when two source ids pool onto one target).
        n_source_keys / n_result_keys: keys carrying ≥1 value before and
            after. The difference is the keys that leave the layer entirely
            because none of their values mapped.
        n_expanded_annotations: input pairs that produced more than one target.
        unmapped_values: the source ids that mapped to nothing, most-used first.
    """

    n_source_values: int = 0
    n_mapped_values: int = 0
    n_source_annotations: int = 0
    n_mapped_annotations: int = 0
    n_result_values: int = 0
    n_result_annotations: int = 0
    n_source_keys: int = 0
    n_result_keys: int = 0
    n_expanded_annotations: int = 0
    unmapped_values: List[str] = field(default_factory=list)

    @property
    def value_coverage(self) -> float:
        """Fraction of distinct source ids that map to a target."""
        if not self.n_source_values:
            return 0.0
        return self.n_mapped_values / self.n_source_values

    @property
    def annotation_coverage(self) -> float:
        """Fraction of key→value annotations that survive the re-keying."""
        if not self.n_source_annotations:
            return 0.0
        return self.n_mapped_annotations / self.n_source_annotations


def remap_values(
    annotations: Dict[str, Set[str]],
    mapping: SupportsTargets,
    label: str,
    *,
    key_label: str = "key",
    value_label: str = "value",
    target_label: str = "target",
) -> Tuple[Dict[str, Set[str]], RemapCoverage]:
    """Re-key the *values* of a ``{key: {source id}}`` map through ``mapping``.

    One-to-many source ids expand to every target; unmapped ones are dropped
    (and counted); keys left with no value at all disappear from the map, as
    they must — the Fisher engine treats an absent protein as having no
    annotation, which is exactly right for a protein whose only ids fall
    outside the target vocabulary.

    ``key_label``/``value_label``/``target_label`` name the axes in the log
    lines (e.g. protein/term/"Disease Ontology term" for the DOID layer).

    Returns:
        ``(remapped map, coverage)``. The coverage report is returned rather
        than only logged so callers and tests can assert on it.
    """
    result: Dict[str, Set[str]] = {}
    value_use: Dict[str, int] = defaultdict(int)
    n_source_annotations = 0
    n_mapped_annotations = 0
    n_expanded = 0

    for key, values in annotations.items():
        mapped: Set[str] = set()
        for value in values:
            n_source_annotations += 1
            value_use[value] += 1
            targets = mapping.targets(value)
            if targets:
                n_mapped_annotations += 1
                if len(targets) > 1:
                    n_expanded += 1
                mapped |= targets
        if mapped:
            result[key] = mapped

    unmapped = sorted(
        (value for value in value_use if not mapping.targets(value)),
        key=lambda value: (-value_use[value], value),
    )
    coverage = RemapCoverage(
        n_source_values=len(value_use),
        n_mapped_values=len(value_use) - len(unmapped),
        n_source_annotations=n_source_annotations,
        n_mapped_annotations=n_mapped_annotations,
        n_result_values=len({value for values in result.values() for value in values}),
        n_result_annotations=sum(len(values) for values in result.values()),
        n_source_keys=len(annotations),
        n_result_keys=len(result),
        n_expanded_annotations=n_expanded,
        unmapped_values=unmapped,
    )

    logger.info(f"Re-keyed annotations {label}:")
    logger.info(
        f"  {value_label} coverage: {coverage.n_mapped_values:,} / "
        f"{coverage.n_source_values:,} distinct source ids "
        f"({coverage.value_coverage:.1%})"
    )
    logger.info(
        f"  annotation coverage: {coverage.n_mapped_annotations:,} / "
        f"{coverage.n_source_annotations:,} {key_label}-{value_label} "
        f"annotations ({coverage.annotation_coverage:.1%})"
    )
    logger.info(
        f"  {key_label}s: {coverage.n_source_keys:,} → "
        f"{coverage.n_result_keys:,}; {value_label}s: "
        f"{coverage.n_source_values:,} → {coverage.n_result_values:,}; "
        f"annotations: {coverage.n_source_annotations:,} → "
        f"{coverage.n_result_annotations:,}"
    )
    if coverage.n_expanded_annotations:
        logger.info(
            f"  one-to-many expansions applied: "
            f"{coverage.n_expanded_annotations:,} annotations"
        )
    if unmapped:
        logger.warning(
            f"  {len(unmapped):,} source ids had no {target_label} and "
            f"were dropped (covering "
            f"{coverage.n_source_annotations - coverage.n_mapped_annotations:,} "
            "annotations); most used: "
            + ", ".join(f"{value} ×{value_use[value]}" for value in unmapped[:5])
        )
    return result, coverage
