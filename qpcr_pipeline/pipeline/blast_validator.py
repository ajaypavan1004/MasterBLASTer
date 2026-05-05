"""
blast_validator.py
Local BLAST+ validator — runs blastn subprocess against a small custom db.

Flow:
  1. build_local_db() — fetches sequences from NCBI, writes FASTA, runs makeblastdb
  2. validate_primer_set() — runs blastn locally, parses XML, checks specificity
  3. validate_all_sets() — pre-filters then BLASTs top N sets

Fallback to remote NCBI BLAST if local db not available.
"""

import logging
import os
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import requests
from Bio import Entrez, SeqIO

from .primer_designer import PrimerSet
from . import config

logger = logging.getLogger(__name__)

# In-process cache: sequence → BlastResult
_blast_cache: Dict[str, "BlastResult"] = {}

# Local db path (set by build_local_db or auto-detected)
_local_db_path: Optional[str] = None


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class BlastHit:
    rank:      int
    accession: str
    title:     str
    organism:  str
    identity:  float
    coverage:  float
    evalue:    float
    bitscore:  float
    is_target: bool = False


@dataclass
class BlastResult:
    query_seq:   str
    query_label: str
    hits:        List[BlastHit] = field(default_factory=list)
    specific:    Optional[bool] = None
    error:       Optional[str]  = None


# ── Pre-filter ────────────────────────────────────────────────────────────────

def _has_homopolymer(seq: str, max_run: int = 4) -> bool:
    seq = seq.upper()
    count = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            count += 1
            if count > max_run:
                return True
        else:
            count = 1
    return False


def prefilter_primer_set(ps: PrimerSet) -> Tuple[bool, List[str]]:
    issues = []
    for name, seq, tm, gc in [
        ("FWD", ps.fwd_seq, ps.fwd_tm, ps.fwd_gc),
        ("REV", ps.rev_seq, ps.rev_tm, ps.rev_gc),
    ]:
        if not (config.PRIMER_MIN_TM <= tm <= config.PRIMER_MAX_TM):
            issues.append(f"{name} Tm {tm:.1f} out of range")
        if not (config.PRIMER_MIN_GC <= gc <= config.PRIMER_MAX_GC):
            issues.append(f"{name} GC {gc:.1f}% out of range")
        if _has_homopolymer(seq):
            issues.append(f"{name} has homopolymer run >4")
    if ps.probe_seq:
        if not (config.PROBE_MIN_GC <= ps.probe_gc <= config.PROBE_MAX_GC):
            issues.append(f"PRB GC {ps.probe_gc:.1f}% out of range")
        if _has_homopolymer(ps.probe_seq):
            issues.append(f"PRB has homopolymer run >4")
    return len(issues) == 0, issues


# ── Local BLAST db builder ────────────────────────────────────────────────────

def build_local_db(
    organism: str,
    db_dir: str = None,
    email: str = "",
    api_key: str = "",
    n_target: int = 10,
    n_relatives: int = 20,
) -> Optional[str]:
    """
    Fetch target + relative sequences from NCBI, build a local BLAST db.
    Returns path to db prefix (for use with blastn -db) or None on failure.
    """
    global _local_db_path

    db_dir = db_dir or os.path.expanduser("~/blast_db")
    os.makedirs(db_dir, exist_ok=True)

    slug = organism.lower().replace(" ", "_")[:30]
    fasta_path = os.path.join(db_dir, f"{slug}_db.fasta")
    db_prefix  = os.path.join(db_dir, f"{slug}_db")

    # If db already exists, reuse it
    if os.path.exists(db_prefix + ".nsi") or os.path.exists(db_prefix + ".nin"):
        logger.info("Local BLAST db already exists: %s", db_prefix)
        _local_db_path = db_prefix
        return db_prefix

    Entrez.email   = email or config.ENTREZ_EMAIL
    Entrez.tool    = config.ENTREZ_TOOL
    if api_key:
        Entrez.api_key = api_key

    from .genome_fetcher import search_taxon_accessions, fetch_sequences, get_close_relatives

    all_records = []

    # Fetch target sequences using the robust genome_fetcher logic
    logger.info("Fetching %d target sequences for local db...", n_target)
    ids = search_taxon_accessions(organism, max_seqs=n_target, email=email, api_key=api_key)
    if ids:
        recs = fetch_sequences(ids, email=email, api_key=api_key)
        all_records.extend(recs)
        logger.info("  Got %d target sequences", len(recs))

    # Fetch relative sequences
    logger.info("Fetching relative sequences for local db...")
    rel_recs = get_close_relatives(organism, n=n_relatives, email=email, api_key=api_key)
    all_records.extend(rel_recs)
    if rel_recs:
        logger.info("  Got %d relative sequences", len(rel_recs))

    if not all_records:
        logger.error("No sequences fetched — cannot build local db")
        return None

    # Write combined FASTA
    with open(fasta_path, "w") as fh:
        SeqIO.write(all_records, fh, "fasta")
    logger.info("Wrote %d sequences to %s", len(all_records), fasta_path)

    # Run makeblastdb
    cmd = [
        "makeblastdb",
        "-in",    fasta_path,
        "-dbtype","nucl",
        "-out",   db_prefix,
        "-title", slug,
    ]
    logger.info("Running makeblastdb...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            logger.error("makeblastdb failed: %s", result.stderr)
            return None
        logger.info("Local BLAST db built: %s", db_prefix)
    except Exception as exc:
        logger.error("makeblastdb error: %s", exc)
        return None

    _local_db_path = db_prefix
    return db_prefix


# ── Local blastn runner ───────────────────────────────────────────────────────

def _run_local_blast(sequences: Dict[str, str], db_path: str) -> Dict[str, List[BlastHit]]:
    """
    Write sequences to a temp FASTA, run blastn locally, parse XML output.
    Returns {label: [BlastHit]} 
    """
    results: Dict[str, List[BlastHit]] = {label: [] for label in sequences}

    # Write query FASTA to temp file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as qf:
        for label, seq in sequences.items():
            qf.write(f">{label}\n{seq}\n")
        query_path = qf.name

    out_path = query_path + ".xml"

    cmd = [
        "blastn",
        "-task",        "blastn-short",
        "-db",          db_path,
        "-query",       query_path,
        "-out",         out_path,
        "-outfmt",      "5",          # XML
        "-word_size",   "7",
        "-evalue",      "1000",
        "-dust",        "no",         # equivalent of FILTER=F for short seqs
        "-num_alignments", "10",
        "-num_threads", "2",
    ]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if proc.returncode != 0:
            logger.error("blastn failed: %s", proc.stderr)
            return results
    except subprocess.TimeoutExpired:
        logger.error("blastn timed out")
        return results
    except FileNotFoundError:
        logger.error("blastn not found — is BLAST+ installed?")
        return results
    finally:
        try:
            os.unlink(query_path)
        except Exception:
            pass

    # Parse XML
    try:
        tree = ET.parse(out_path)
        root = tree.getroot()
    except Exception as exc:
        logger.error("XML parse error: %s", exc)
        return results
    finally:
        try:
            os.unlink(out_path)
        except Exception:
            pass

    labels = list(sequences.keys())
    for iteration in root.findall(".//Iteration"):
        query_def = iteration.findtext("Iteration_query-def", "").strip()
        query_len = int(iteration.findtext("Iteration_query-len", "1"))

        label = None
        for l in labels:
            if l in query_def or query_def in l:
                label = l
                break
        if label is None:
            label = query_def

        hits: List[BlastHit] = []
        for rank, hit_el in enumerate(iteration.findall(".//Hit"), start=1):
            title     = hit_el.findtext("Hit_def", "")
            accession = hit_el.findtext("Hit_accession", "")
            hsp       = hit_el.find(".//Hsp")
            if hsp is None:
                continue
            align_len    = int(hsp.findtext("Hsp_align-len", "1"))
            identity     = int(hsp.findtext("Hsp_identity",  "0"))
            evalue       = float(hsp.findtext("Hsp_evalue",   "999"))
            bitscore     = float(hsp.findtext("Hsp_bit-score","0"))
            pct_identity = (identity / align_len * 100) if align_len else 0
            coverage     = (align_len / query_len * 100) if query_len else 0
            organism     = " ".join(title.split()[:2]) if title else "Unknown"
            hits.append(BlastHit(
                rank=rank, accession=accession, title=title,
                organism=organism, identity=pct_identity,
                coverage=coverage, evalue=evalue, bitscore=bitscore,
            ))
            if rank >= 10:
                break

        results[label] = hits

    return results


def _evaluate(hits: List[BlastHit], target_organism: str) -> Tuple[bool, List[BlastHit]]:
    target_tokens = set(target_organism.lower().split())
    for hit in hits:
        hit.is_target = bool(target_tokens & set(hit.organism.lower().split()))
    return all(h.is_target for h in hits), hits


# ── Public API ────────────────────────────────────────────────────────────────

def validate_primer_set(
    ps: PrimerSet,
    target_organism: str,
    db_path: str = None,
) -> PrimerSet:
    """
    BLAST fwd + rev + probe in one shot using local BLAST+.
    Uses cache to skip re-blasting identical sequences.
    """
    db = db_path or _local_db_path
    if not db:
        logger.error("No local BLAST db available. Run build_local_db() first.")
        ps.blast_pass = None
        return ps

    to_blast: Dict[str, str] = {}
    cached:   Dict[str, BlastResult] = {}

    for label, seq in [
        ("forward_primer", ps.fwd_seq),
        ("reverse_primer", ps.rev_seq),
        ("probe",          ps.probe_seq),
    ]:
        if not seq:
            continue
        if seq in _blast_cache:
            logger.info("Cache hit for %s", label)
            cached[label] = _blast_cache[seq]
        else:
            to_blast[label] = seq

    if to_blast:
        logger.info("Running local BLAST for %s", list(to_blast.keys()))
        hit_map = _run_local_blast(to_blast, db)

        for label, seq in to_blast.items():
            hits = hit_map.get(label, [])
            specific, hits = _evaluate(hits, target_organism)
            r = BlastResult(query_seq=seq, query_label=label, hits=hits, specific=specific)
            _blast_cache[seq] = r
            cached[label] = r
            logger.info("  %s: %d hits, specific=%s", label, len(hits), specific)

    ps.fwd_blast   = cached.get("forward_primer")
    ps.rev_blast   = cached.get("reverse_primer")
    ps.probe_blast = cached.get("probe")

    checks = [r for r in [ps.fwd_blast, ps.rev_blast, ps.probe_blast]
              if r is not None and r.error is None]
    ps.blast_pass = all(r.specific for r in checks) if checks else None
    return ps


def validate_all_sets(
    primer_sets: List[PrimerSet],
    target_organism: str,
    max_sets: int = 3,
    db_path: str = None,
) -> List[PrimerSet]:
    db = db_path or _local_db_path
    if not db:
        logger.error("No local BLAST db. Call build_local_db() before validate_all_sets().")
        return primer_sets

    # Pre-filter
    filtered = []
    for ps in primer_sets:
        passed, issues = prefilter_primer_set(ps)
        if passed:
            filtered.append(ps)
        else:
            logger.debug("Pre-filter rejected pair #%d: %s", ps.pair_index, issues)

    if not filtered:
        logger.warning("All sets failed pre-filter — blasting top %d anyway", max_sets)
        filtered = primer_sets

    to_blast = filtered[:max_sets]
    logger.info("Local BLAST: %d set(s) against '%s'", len(to_blast), target_organism)

    for i, ps in enumerate(to_blast):
        logger.info("--- BLAST run %d/%d (Pair #%d) ---", i + 1, len(to_blast), ps.pair_index)
        validate_primer_set(ps, target_organism, db_path=db)

    return primer_sets
