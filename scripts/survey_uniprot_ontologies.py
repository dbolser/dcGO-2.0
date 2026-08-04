#!/usr/bin/env python
"""Survey which vocabularies the UniProt flat file can supply as ontologies.

Adding an ontology to dcGO is cheap when UniProt already carries it per
accession (no identifier mapping needed) — but only if it *is* a vocabulary.
Most of the ~150 ``DR`` databases in Swiss-Prot are 1:1 mirrors of the protein
(AlphaFoldDB, STRING, GeneCards …); as "terms" they would give every
contingency table a single protein and no signal at all.

This script measures the difference: for one taxon, how many proteins and how
many *distinct* terms each ``DR`` database (and the ``KW`` keyword vocabulary)
contributes. A high proteins-per-term ratio means a usable ontology; ~1.0 means
an accession mirror. The result is the evidence behind the registry in
``src/ontology_registry.py`` — see ``docs/uniprot_ontology_survey.md``.

Usage
-----
    uv run python scripts/survey_uniprot_ontologies.py --output docs/dr_survey.tsv
    uv run python scripts/survey_uniprot_ontologies.py --taxid 10090   # mouse
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.uniprot_annotation_source import _open_text  # noqa: E402

#: ``OX   NCBI_TaxID=9606 {ECO:0000313|EMBL:ABC12345.1};`` — the evidence tag is
#: stripped before the taxon id is compared.
_OX_EVIDENCE = re.compile(r"\{[^}]*\}")

#: A ``DR`` third field is only a *type* if the database uses a handful of
#: values (``MIM``'s gene/phenotype, ``Pharos``'s four development levels). Most
#: databases put another identifier or a free-text name there, which would
#: fragment the report into one row per protein — so the typed breakdown is only
#: emitted for databases below this many distinct values.
MAX_ID_TYPES_PER_DATABASE = 20


def survey(path: Path, taxid: str) -> Tuple[Dict, Dict, Dict, int, int]:
    """Count proteins and distinct terms per ``(database, id_type)`` for one taxon.

    Returns ``(proteins, terms, examples, n_entries, n_taxon_entries)`` where the
    first three are keyed by ``(database, id_type)``; the empty ``id_type`` key
    holds the database total.
    """
    proteins: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    terms: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    examples: Dict[Tuple[str, str], str] = {}
    n_entries = n_taxon = 0

    accession: str | None = None
    in_taxon = False
    cross_refs = []
    keyword_lines: list[str] = []
    marker = f"NCBI_TaxID={taxid}"

    with _open_text(path) as handle:
        for line in handle:
            tag = line[:2]
            if tag == "AC":
                if accession is None:
                    accession = line[5:].split(";")[0].strip() or None
            elif tag == "OX":
                # Compare whole semicolon-delimited fields, minus any
                # {ECO:…} evidence tag: a substring test would let
                # NCBI_TaxID=96060 match a query for taxon 9606.
                in_taxon = in_taxon or any(
                    _OX_EVIDENCE.sub("", field).strip().rstrip(".") == marker
                    for field in line[5:].split(";")
                )
            elif tag == "DR":
                fields = line[5:].split(";")
                if len(fields) >= 2:
                    database = fields[0].strip()
                    external_id = fields[1].strip()
                    id_type = fields[2].strip().rstrip(".") if len(fields) >= 3 else ""
                    if database and external_id:
                        cross_refs.append((database, external_id, id_type))
            elif tag == "KW":
                keyword_lines.append(line[5:])
            elif line.startswith("//"):
                n_entries += 1
                if in_taxon and accession:
                    n_taxon += 1
                    for database, external_id, id_type in cross_refs:
                        # Database total, then the typed breakdown; the latter is
                        # dropped at output time for databases that turn out to
                        # have too many distinct third-field values to be types.
                        for key in ((database, ""), (database, id_type)):
                            proteins[key].add(accession)
                            terms[key].add(external_id)
                            examples.setdefault(key, external_id)
                    text = " ".join(p.strip() for p in keyword_lines).rstrip(" .")
                    for keyword in (k.strip() for k in text.split(";")):
                        if keyword:
                            proteins[("KW", "")].add(accession)
                            terms[("KW", "")].add(keyword)
                            examples.setdefault(("KW", ""), keyword)
                accession, in_taxon = None, False
                cross_refs, keyword_lines = [], []

    return proteins, terms, examples, n_entries, n_taxon


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--uniprot-dat",
        type=Path,
        default=Path("data/raw/uniprot_sprot_dat/uniprot_sprot.dat.gz"),
    )
    parser.add_argument(
        "--taxid", default="9606", help="NCBI taxon id (default: human)"
    )
    parser.add_argument("--output", type=Path, default=None, help="TSV output path")
    args = parser.parse_args()

    proteins, terms, examples, n_entries, n_taxon = survey(args.uniprot_dat, args.taxid)

    # Keep the typed rows only where the third field behaves like a type.
    id_types_per_database: Dict[str, int] = defaultdict(int)
    for database, id_type in proteins:
        if id_type:
            id_types_per_database[database] += 1
    reportable = {
        key
        for key in proteins
        if not key[1] or id_types_per_database[key[0]] <= MAX_ID_TYPES_PER_DATABASE
    }

    lines = [
        f"# entries={n_entries} taxon={args.taxid} taxon_entries={n_taxon}",
        "database\tid_type\tn_proteins\tn_distinct_terms\tproteins_per_term\texample",
    ]
    for key in sorted(reportable, key=lambda k: (-len(proteins[k]), k[0], k[1])):
        n_proteins, n_terms = len(proteins[key]), len(terms[key])
        ratio = n_proteins / n_terms if n_terms else 0.0
        lines.append(
            f"{key[0]}\t{key[1]}\t{n_proteins}\t{n_terms}\t{ratio:.2f}\t{examples[key]}"
        )

    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
        print(f"Wrote {len(lines) - 2} rows to {args.output}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
