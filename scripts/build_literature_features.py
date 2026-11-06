"""Per-protein literature embedding features from BioEmbed (MedEmbed-small).

Each human protein gets one 384-d vector summarising the PubMed literature that
mentions its gene: the weighted mean of the MedEmbed-small title+abstract
vectors of every paper PubTator3 tagged with that gene, L2-normalised.

**Temporal rule.** This feeds a temporal benchmark, so for a cutoff year ``Y``
only papers whose publication year is ``<= Y`` contribute — a paper published
after the cutoff is never read, at any stage (``n_papers``/``w_sum`` count the
same papers as ``emb``). Papers whose year is unknown (pmid absent from the
abstracts index, or a null year) are excluded at every cutoff and counted in
the sidecar JSON, never silently assumed old.

**Weighting.** A paper tagged with ``n`` human genes gives weight ``1/n`` to
each of them, so a single-gene paper counts fully while a 10,000-gene screen
barely moves any one gene. ``n`` is the paper's full HGNC gene count, including
genes that map to no protein, so a weight does not depend on the mapping.

**Gene → protein.** NCBI GeneID → UniProt accession, from two sources unioned:
the Swiss-Prot flat file's ``DR GeneID`` lines
(:func:`src.gene_mapping.parse_gene_accession_index`, restricted to
``OX NCBI_TaxID=9606`` entries) and ``idmapping_selected.tab`` rows whose
NCBI-taxon column is 9606, which add the human TrEMBL accessions. A gene's
papers go to every accession mapped to it (one-to-many is kept, as in
:mod:`src.gene_mapping`). An accession mapped to several genes (readthrough
loci) takes the *union* of their papers, each paper counted once.

Paper vectors are stored as float16, so their norms drift by up to ~0.5%; each
is re-normalised in float32 before averaging so every paper carries exactly its
weight. Accessions with no paper ``<= Y`` are omitted from that cutoff's file.

Outputs, per cutoff, in ``--output-dir``:

* ``medembed_protein_y<Y>.npz`` — ``accessions`` (str), ``emb`` (float32,
  n × 384), ``n_papers`` (int32), ``w_sum`` (float32).
* ``medembed_protein_y<Y>.json`` — inputs with SHA-256, counts, mapping
  ambiguity and GOA coverage.

pyarrow is not a project dependency, so run with::

    uv run --with pyarrow python scripts/build_literature_features.py \\
        --bioembed-root /path/to/BioEmbed --cutoff-year 2017 --cutoff-year 2020
"""

from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from multiprocessing.pool import Pool
from pathlib import Path
from typing import BinaryIO, Dict, Iterator, List, Set, Tuple

import numpy as np
import pyarrow.compute as pc
import pyarrow.parquet as pq
import scipy.sparse as sp
from loguru import logger

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.gene_mapping import parse_gene_accession_index  # noqa: E402
from src.hierarchy import open_text  # noqa: E402
from src.run_manifest import sha256_file, utc_now  # noqa: E402

HUMAN_TAXON = "9606"
DEFAULT_CUTOFFS = [2017, 2020]
DEFAULT_GAFS = [
    REPO / "data/raw/goa_archive/goa_human.gaf.181.gz",
    REPO / "data/raw/goa_archive/goa_human.gaf.205.gz",
    REPO / "data/raw/goa_annotations/goa_human.gaf.gz",
]
# Layout inside --bioembed-root.
GENE_PMIDS = "data/geometry/human_gene_pmids.parquet"
VECTORS = "data/geometry/medembed_small/X.npy"
VECTOR_PMIDS = "data/geometry/medembed_small/pmids.npy"
ABSTRACTS = "genesearch/index/abstracts.parquet"
# idmapping_selected.tab columns (0-based): 0 UniProtKB-AC, 2 GeneID, 12 NCBI-taxon.
IDMAP_ACC, IDMAP_GENEID, IDMAP_TAXON = 0, 2, 12
_OX = re.compile(r"OX   NCBI_TaxID=(\d+)")


# --------------------------------------------------------------------------- #
# Gene → accession maps (run in worker processes)
# --------------------------------------------------------------------------- #


def taxon_accessions(dat_path: Path, taxon: str) -> Set[str]:
    """Primary accessions of the flat-file entries whose ``OX`` taxon is ``taxon``."""
    accessions: Set[str] = set()
    accession: str | None = None
    with open_text(dat_path, label="UniProt flat file") as handle:
        for line in handle:
            if line.startswith("AC   ") and accession is None:
                accession = line[5:].split(";")[0].strip() or None
            elif line.startswith("OX   ") and accession is not None:
                match = _OX.match(line)
                if match and match.group(1) == taxon:
                    accessions.add(accession)
            elif line.startswith("//"):
                accession = None
    return accessions


def sprot_geneid_map(dat_path: Path) -> Dict[str, Set[str]]:
    """Swiss-Prot ``DR GeneID`` map, all species (taxon filter applied later)."""
    return parse_gene_accession_index(dat_path).geneid.source_to_accessions


@contextmanager
def open_binary(path: Path) -> Iterator[BinaryIO]:
    """Binary read handle; ``.gz`` is piped through ``pigz`` when installed."""
    pigz = shutil.which("pigz")
    if path.suffix != ".gz":
        with path.open("rb") as handle:
            yield handle
    elif pigz is None:
        with gzip.open(path, "rb") as handle:
            yield handle  # type: ignore[misc]
    else:
        proc = subprocess.Popen([pigz, "-dc", str(path)], stdout=subprocess.PIPE)
        assert proc.stdout is not None
        try:
            yield proc.stdout
        except BaseException:
            proc.kill()
            proc.wait()
            raise
        proc.stdout.close()
        if proc.wait() != 0:
            raise RuntimeError(f"pigz failed on {path} (exit {proc.returncode})")


def idmapping_geneid_map(
    path: Path, taxon: str
) -> Tuple[Dict[str, Set[str]], Set[str]]:
    """GeneID → accessions from ``idmapping_selected.tab`` rows of one taxon.

    Also returns every accession of that taxon, GeneID or not — the current
    UniProtKB proteome the coverage figures are checked against.

    The file is ~30 GB decompressed and ``taxon`` rows are a tiny fraction, so
    instead of splitting every line this searches 64 MB blocks for
    ``\\t<taxon>\\t`` and only splits the lines that hit, then confirms the hit
    is in the NCBI-taxon column. The GeneID column holds ``"; "``-separated ids.
    """
    want = taxon.encode()
    needle = b"\t" + want + b"\t"
    mapping: Dict[str, Set[str]] = defaultdict(set)
    accessions: Set[str] = set()

    def scan(buf: bytes) -> None:
        start = 0
        while (hit := buf.find(needle, start)) != -1:
            line_start = buf.rfind(b"\n", 0, hit) + 1
            line_end = buf.find(b"\n", hit)
            line_end = len(buf) if line_end == -1 else line_end
            fields = buf[line_start:line_end].split(b"\t")
            if len(fields) > IDMAP_TAXON and fields[IDMAP_TAXON] == want:
                accession = fields[IDMAP_ACC].decode()
                accessions.add(accession)
                for gene in fields[IDMAP_GENEID].decode().split(";"):
                    if gene.strip():
                        mapping[gene.strip()].add(accession)
            start = line_end + 1

    tail = b""
    with open_binary(path) as handle:
        while block := handle.read(1 << 26):
            buf = tail + block
            cut = buf.rfind(b"\n") + 1
            scan(buf[:cut])
            tail = buf[cut:]
    scan(tail)
    return dict(mapping), accessions


def gaf_accessions(path: Path) -> Set[str]:
    """UniProtKB accessions (column 2) of a GAF file."""
    accessions: Set[str] = set()
    with open_text(path, label="GAF file") as handle:
        for line in handle:
            if not line.startswith("!"):
                db, accession = line.split("\t", 2)[:2]
                if db == "UniProtKB":
                    accessions.add(accession)
    return accessions


# --------------------------------------------------------------------------- #
# Papers
# --------------------------------------------------------------------------- #


def lookup(keys: np.ndarray, table_keys: np.ndarray) -> np.ndarray:
    """Row index of each of ``keys`` in unsorted unique ``table_keys``, -1 if absent."""
    order = np.argsort(table_keys, kind="stable")
    sorted_keys = table_keys[order]
    pos = np.searchsorted(sorted_keys, keys).clip(max=len(sorted_keys) - 1)
    return np.where(sorted_keys[pos] == keys, order[pos], -1)


def load_vectors(path: Path) -> np.ndarray:
    """float16 paper vectors → re-normalised float32, converted in row blocks."""
    raw = np.load(path, mmap_mode="r")
    out = np.empty(raw.shape, dtype=np.float32)
    step = 1 << 20

    def convert(lo: int) -> None:
        block = np.asarray(raw[lo : lo + step], dtype=np.float32)
        block /= np.linalg.norm(block, axis=1, keepdims=True)
        out[lo : lo + step] = block

    with ThreadPoolExecutor(max_workers=16) as pool:
        list(pool.map(convert, range(0, len(raw), step)))
    return out


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--cutoff-year",
        type=int,
        action="append",
        help=f"keep papers with year <= this; repeatable (default {DEFAULT_CUTOFFS})",
    )
    parser.add_argument(
        "--bioembed-root",
        type=Path,
        required=True,
        help=f"BioEmbed checkout holding {GENE_PMIDS}, {VECTORS} and {ABSTRACTS}",
    )
    parser.add_argument(
        "--sprot",
        type=Path,
        default=REPO / "data/raw/uniprot_sprot_dat/uniprot_sprot.dat.gz",
    )
    parser.add_argument(
        "--idmapping",
        type=Path,
        default=REPO / "data/raw/uniprot_idmapping/idmapping_selected.tab.gz",
    )
    parser.add_argument(
        "--coverage-gaf",
        type=Path,
        action="append",
        help="GAF whose accessions coverage is reported for; repeatable "
        "(default: the 2018 and 2021 t0 and the current human GOA)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=REPO / "data/processed/literature"
    )
    args = parser.parse_args(argv)
    args.cutoff_year = sorted(set(args.cutoff_year or DEFAULT_CUTOFFS))
    args.coverage_gaf = args.coverage_gaf or DEFAULT_GAFS
    return args


def main(argv: List[str] | None = None) -> None:
    args = parse_args(argv)
    started = time.monotonic()
    root = args.bioembed_root
    inputs = {
        "gene_pmids": root / GENE_PMIDS,
        "vectors": root / VECTORS,
        "vector_pmids": root / VECTOR_PMIDS,
        "abstracts": root / ABSTRACTS,
        "sprot": args.sprot,
        "idmapping": args.idmapping,
        **{f"gaf:{path.name}": path for path in args.coverage_gaf},
    }
    for name, path in inputs.items():
        if not path.exists():
            raise FileNotFoundError(f"{name} input not found: {path}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # The three UniProt scans are pure-Python and independent: run them in
    # worker processes while this process hashes inputs and loads papers.
    # Pool's context manager terminates the workers, so a failure exits at
    # once instead of waiting for the scans to finish.
    with (
        multiprocessing.get_context("spawn").Pool(3) as pool,
        ThreadPoolExecutor(len(inputs)) as hashers,
    ):
        build(args, inputs, pool, hashers)
    logger.info(f"Done in {time.monotonic() - started:.0f} s")


def build(
    args: argparse.Namespace,
    inputs: Dict[str, Path],
    pool: Pool,
    hashers: ThreadPoolExecutor,
) -> None:
    """Everything after argument checks: map genes, then write each cutoff."""
    created = utc_now()
    r_human = pool.apply_async(taxon_accessions, (args.sprot, HUMAN_TAXON))
    r_sprot = pool.apply_async(sprot_geneid_map, (args.sprot,))
    r_idmap = pool.apply_async(idmapping_geneid_map, (args.idmapping, HUMAN_TAXON))
    f_hashes = {name: hashers.submit(sha256_file, p) for name, p in inputs.items()}

    # ---- papers: genes, weights, vector rows, years ----------------------- #
    logger.info(f"Reading {inputs['gene_pmids']}")
    table = pq.read_table(inputs["gene_pmids"])
    pmid = table["pmid"].to_numpy()
    n_human_genes = table["n_human_genes"].to_numpy()
    gene_lists = table["ncbi_genes"].combine_chunks()
    list_len = pc.list_value_length(gene_lists).to_numpy(zero_copy_only=False)
    if len(np.unique(pmid)) != len(pmid):
        raise ValueError("duplicate pmids in the gene table")
    if not np.array_equal(list_len, n_human_genes):
        raise ValueError("n_human_genes disagrees with len(ncbi_genes)")
    encoded = pc.dictionary_encode(pc.list_flatten(gene_lists))
    gene_names: List[str] = encoded.dictionary.to_pylist()
    pair_gene = encoded.indices.to_numpy().astype(np.int64)
    pair_paper = pc.list_parent_indices(gene_lists).to_numpy().astype(np.int64)
    n_dup_pairs = len(pair_gene) - len(
        np.unique(pair_paper * len(gene_names) + pair_gene)
    )
    if n_dup_pairs:
        raise ValueError(f"{n_dup_pairs} duplicate (paper, gene) pairs")
    logger.info(
        f"  {len(pmid):,} papers, {len(pair_gene):,} paper-gene pairs, "
        f"{len(gene_names):,} genes"
    )

    vector_pmids = np.load(inputs["vector_pmids"])
    vector_rows = np.load(inputs["vectors"], mmap_mode="r").shape[0]
    if vector_pmids.shape != (vector_rows,):
        raise ValueError(f"pmids.npy {vector_pmids.shape} vs X.npy rows {vector_rows}")
    if len(np.unique(vector_pmids)) != len(vector_pmids):
        raise ValueError("duplicate pmids in pmids.npy")
    vec_row = lookup(pmid, vector_pmids)

    logger.info(f"Reading pmid, year from {inputs['abstracts']}")
    abstracts = pq.read_table(inputs["abstracts"], columns=["pmid", "year"])
    abs_pmid = abstracts["pmid"].to_numpy()
    abs_year = pc.fill_null(abstracts["year"], -1).to_numpy().astype(np.int32)
    if len(np.unique(abs_pmid)) != len(abs_pmid):
        raise ValueError("duplicate pmids in abstracts.parquet")
    abs_row = lookup(pmid, abs_pmid)
    year = np.where(abs_row >= 0, abs_year[abs_row.clip(min=0)], -1)
    has_vec = vec_row >= 0
    has_year = year >= 0
    paper_counts = {
        "papers_in_gene_table": len(pmid),
        "paper_gene_pairs": len(pair_gene),
        "genes_in_gene_table": len(gene_names),
        "papers_without_vector": int((~has_vec).sum()),
        "papers_with_vector_not_in_abstracts": int((has_vec & (abs_row < 0)).sum()),
        "papers_with_vector_null_year": int(
            (has_vec & (abs_row >= 0) & ~has_year).sum()
        ),
        "papers_with_vector_without_year": int((has_vec & ~has_year).sum()),
        "year_range_of_usable_papers": [
            int(year[has_vec & has_year].min()),
            int(year[has_vec & has_year].max()),
        ],
    }
    logger.info(f"  {paper_counts}")
    vectors = load_vectors(inputs["vectors"])
    logger.info(f"  vectors {vectors.shape} loaded")

    # ---- gene → accession ------------------------------------------------- #
    human: Set[str] = r_human.get()
    sprot_all: Dict[str, Set[str]] = r_sprot.get()
    idmap, idmap_human = r_idmap.get()
    genes_here = set(gene_names)
    sprot_map: Dict[str, Set[str]] = {}
    n_nonhuman_dropped = 0
    for gene in genes_here & sprot_all.keys():
        kept = sprot_all[gene] & human
        n_nonhuman_dropped += len(sprot_all[gene]) - len(kept)
        if kept:
            sprot_map[gene] = kept
    gene_to_acc: Dict[str, Set[str]] = defaultdict(set)
    for source in (sprot_map, idmap):
        for gene in genes_here & source.keys():
            gene_to_acc[gene] |= source[gene]
    acc_to_genes: Dict[str, Set[str]] = defaultdict(set)
    for gene, accs in gene_to_acc.items():
        for acc in accs:
            acc_to_genes[acc].add(gene)
    sprot_accs = set().union(*sprot_map.values())
    all_accs = sorted(acc_to_genes)
    fan_out = np.array([len(a) for a in gene_to_acc.values()])
    pairs_mapped = np.isin(
        pair_gene, [i for i, g in enumerate(gene_names) if g in gene_to_acc]
    )
    mapping = {
        "genes_mapped_swissprot": len(sprot_map),
        "genes_mapped_idmapping": len(genes_here & idmap.keys()),
        "genes_mapped_union": len(gene_to_acc),
        "genes_mapped_only_via_idmapping": len(gene_to_acc.keys() - sprot_map.keys()),
        "genes_unmapped": len(genes_here) - len(gene_to_acc),
        "paper_gene_pairs_on_unmapped_genes": int((~pairs_mapped).sum()),
        "swissprot_nonhuman_accessions_dropped": n_nonhuman_dropped,
        "accessions_mapped": len(all_accs),
        "accessions_reviewed": sum(a in human for a in all_accs),
        "accessions_unreviewed": sum(a not in human for a in all_accs),
        "accessions_via_swissprot_geneid": len(sprot_accs),
        "genes_one_to_many_swissprot": sum(len(a) > 1 for a in sprot_map.values()),
        "genes_one_to_many_union": int((fan_out > 1).sum()),
        "max_accessions_per_gene": int(fan_out.max()),
        "accessions_many_to_one_gene": sum(len(g) > 1 for g in acc_to_genes.values()),
        "max_genes_per_accession": max(len(g) for g in acc_to_genes.values()),
    }
    logger.info(f"  mapping {mapping}")

    # Accessions sharing a gene set get identical features: compute once per set.
    gene_index = {g: i for i, g in enumerate(gene_names)}
    set_ids: Dict[frozenset[str], int] = {}
    acc_set = np.array(
        [set_ids.setdefault(frozenset(acc_to_genes[a]), len(set_ids)) for a in all_accs]
    )
    rows, cols = [], []
    for gene_set, sid in set_ids.items():
        for gene in gene_set:
            rows.append(sid)
            cols.append(gene_index[gene])
    set_by_gene = sp.csr_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, cols)),
        shape=(len(set_ids), len(gene_names)),
    )
    # Paper weight indexed by vector row, the column space of the incidence.
    weight = np.zeros(len(vectors), dtype=np.float32)
    weight[vec_row[has_vec]] = 1.0 / n_human_genes[has_vec]
    gaf_sets = {path.name: gaf_accessions(path) for path in args.coverage_gaf}
    current_human = human | idmap_human
    hashes = {name: f.result() for name, f in f_hashes.items()}

    # ---- per cutoff ------------------------------------------------------- #
    for cutoff in args.cutoff_year:
        t0 = time.monotonic()
        usable = has_vec & has_year & (year <= cutoff)
        keep = usable[pair_paper]
        gene_by_paper = sp.csr_matrix(
            (
                np.ones(int(keep.sum()), dtype=np.float32),
                (pair_gene[keep], vec_row[pair_paper[keep]]),
            ),
            shape=(len(gene_names), len(vectors)),
        )
        incidence = (set_by_gene @ gene_by_paper).tocsr()
        incidence.sum_duplicates()
        incidence.data = weight[incidence.indices]
        set_papers = np.diff(incidence.indptr).astype(np.int32)
        set_wsum = np.asarray(incidence.sum(axis=1)).ravel().astype(np.float32)
        set_emb = np.asarray(incidence @ vectors, dtype=np.float32)

        present = set_papers[acc_set] > 0
        accs = np.array(all_accs)[present]
        emb = set_emb[acc_set[present]]
        norms = np.linalg.norm(emb, axis=1, keepdims=True)
        if (norms == 0).any():
            raise ValueError("zero-norm literature vector")
        emb /= norms
        n_papers = set_papers[acc_set[present]]
        w_sum = set_wsum[acc_set[present]]

        stem = f"medembed_protein_y{cutoff}"
        npz_path = args.output_dir / f"{stem}.npz"
        tmp = npz_path.with_suffix(".tmp.npz")
        np.savez(tmp, accessions=accs, emb=emb, n_papers=n_papers, w_sum=w_sum)
        tmp.replace(npz_path)

        covered = set(accs.tolist())
        coverage = {}
        for name, gaf in gaf_sets.items():
            coverage[name] = {
                "accessions": len(gaf),
                "in_current_human_uniprot": len(gaf & current_human),
                "mapped_to_a_gene": len(gaf & acc_to_genes.keys()),
                "with_papers": len(gaf & covered),
                "with_papers_reviewed": len(gaf & covered & human),
                "with_papers_only_via_idmapping": len((gaf & covered) - sprot_accs),
                "fraction_with_papers": round(len(gaf & covered) / len(gaf), 4),
            }
        genes_used = np.unique(pair_gene[keep])
        counts = {
            **paper_counts,
            "papers_le_cutoff": int(usable.sum()),
            "papers_le_cutoff_with_mapped_gene": int(
                np.unique(pair_paper[keep & pairs_mapped]).size
            ),
            "papers_after_cutoff_excluded": int((has_vec & has_year & ~usable).sum()),
            "genes_le_cutoff": int(genes_used.size),
            "genes_le_cutoff_mapped": int(
                sum(gene_names[i] in gene_to_acc for i in genes_used)
            ),
            "accessions": int(len(accs)),
            "accessions_mapped_without_papers_le_cutoff": int((~present).sum()),
            "n_papers_median": float(np.median(n_papers)),
            "n_papers_max": int(n_papers.max()),
        }
        sidecar = {
            "created": created,
            "script": "scripts/build_literature_features.py",
            "command": sys.argv,
            "cutoff_year": cutoff,
            "temporal_rule": "papers with publication year <= cutoff_year only; "
            "papers without a year are excluded",
            "weighting": "paper weight 1/n_human_genes; emb = L2-normalised "
            "weighted mean of re-normalised float32 paper vectors; an accession "
            "takes the union of its genes' papers",
            "embedding_model": "MedEmbed-small (BioEmbed), title+abstract, 384-d",
            "inputs": {
                name: {"path": str(path), "sha256": hashes[name]}
                for name, path in inputs.items()
            },
            "output": {"path": str(npz_path), "sha256": sha256_file(npz_path)},
            "counts": counts,
            "mapping": mapping,
            "coverage": coverage,
            "seconds_this_cutoff": round(time.monotonic() - t0, 1),
        }
        json_path = npz_path.with_suffix(".json")
        json_path.write_text(json.dumps(sidecar, indent=2) + "\n")
        logger.info(f"Y={cutoff}: {len(accs):,} accessions -> {npz_path}")
        logger.info(f"  counts {counts}")
        logger.info(f"  coverage {coverage}")


if __name__ == "__main__":
    main()
