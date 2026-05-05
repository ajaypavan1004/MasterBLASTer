"""
blast_validator.py
Submits primer/probe sequences to NCBI BLAST (blastn) via the URL API,
runs N workers in parallel with asyncio, and parses hit organism names
to confirm specificity.
"""

import asyncio
import logging
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import aiohttp
import requests

from .primer_designer import PrimerSet
from . import config

logger = logging.getLogger(__name__)

BLAST_PUT_URL  = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
BLAST_GET_URL  = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class BlastHit:
    rank:        int
    accession:   str
    title:       str
    organism:    str          # extracted from hit title
    identity:    float        # percent identity
    coverage:    float        # query coverage %
    evalue:      float
    bitscore:    float
    is_target:   bool = False # set during evaluation

@dataclass
class BlastResult:
    query_seq:   str
    query_label: str          # "forward_primer" | "reverse_primer" | "probe"
    hits:        List[BlastHit] = field(default_factory=list)
    specific:    Optional[bool] = None   # True if top hits match target
    error:       Optional[str]  = None


# ── NCBI BLAST URL API helpers ────────────────────────────────────────────────

def _submit_blast(sequence: str, session: requests.Session) -> Optional[str]:
    """POST a sequence to BLAST; return RID (request ID) or None."""
    params = {
        "CMD":      "Put",
        "PROGRAM":  config.BLAST_PROGRAM,
        "DATABASE": config.BLAST_DB,
        "QUERY":    sequence,
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": config.BLAST_HITLIST,
        "WORD_SIZE":    config.BLAST_WORD_SIZE,
        "EXPECT":       config.BLAST_EVALUE,
        "FILTER":       "L",          # low complexity filter
        "SHORT_QUERY_ADJUST": "true",
    }
    try:
        resp = session.post(BLAST_PUT_URL, data=params, timeout=30)
        resp.raise_for_status()
    except Exception as exc:
        logger.error("BLAST submit failed: %s", exc)
        return None

    # Extract RID from response
    for line in resp.text.splitlines():
        if line.startswith("    RID = "):
            return line.split("=")[1].strip()
    logger.error("Could not parse RID from BLAST response")
    return None


def _poll_blast(rid: str, session: requests.Session, timeout: int = 120) -> Optional[str]:
    """Poll until BLAST job is READY; return XML string or None."""
    params = {
        "CMD":         "Get",
        "RID":         rid,
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": config.BLAST_HITLIST,
    }
    deadline = time.time() + timeout
    while time.time() < deadline:
        time.sleep(5)
        try:
            resp = session.get(BLAST_GET_URL, params=params, timeout=30)
            resp.raise_for_status()
        except Exception as exc:
            logger.warning("BLAST poll error: %s", exc)
            continue

        status_line = ""
        for line in resp.text.splitlines():
            if "Status=" in line:
                status_line = line.strip()
                break

        if "READY" in status_line:
            return resp.text
        elif "FAILED" in status_line:
            logger.error("BLAST job %s FAILED", rid)
            return None
        # else WAITING — continue polling

    logger.error("BLAST job %s timed out after %ds", rid, timeout)
    return None


def _parse_blast_xml(xml_text: str) -> List[BlastHit]:
    """Parse BLAST XML output into BlastHit objects."""
    hits: List[BlastHit] = []
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        logger.error("XML parse error: %s", exc)
        return hits

    # Navigate: BlastOutput > BlastOutput_iterations > Iteration > Iteration_hits > Hit
    for rank, hit_el in enumerate(root.iter("Hit"), start=1):
        title = hit_el.findtext("Hit_def", "")
        accession = hit_el.findtext("Hit_accession", "")

        # Best HSP
        hsp = hit_el.find(".//Hsp")
        if hsp is None:
            continue

        align_len    = int(hsp.findtext("Hsp_align-len", "1"))
        identity     = int(hsp.findtext("Hsp_identity",  "0"))
        query_len    = int(hit_el.findtext("../Iteration_query-len",
                           hsp.findtext("Hsp_query-to", "1")))
        evalue       = float(hsp.findtext("Hsp_evalue",   "999"))
        bitscore     = float(hsp.findtext("Hsp_bit-score","0"))

        pct_identity = (identity / align_len * 100) if align_len else 0
        coverage     = (align_len / query_len * 100) if query_len else 0

        # Extract organism from title (first word pair after ">" or start)
        # Titles look like: "Cyclospora cayetanensis isolate XY complete genome"
        organism = " ".join(title.split()[:2]) if title else "Unknown"

        hits.append(BlastHit(
            rank=rank,
            accession=accession,
            title=title,
            organism=organism,
            identity=pct_identity,
            coverage=coverage,
            evalue=evalue,
            bitscore=bitscore,
        ))

        if rank >= config.BLAST_HITLIST:
            break

    return hits


def _evaluate_specificity(hits: List[BlastHit], target_organism: str) -> Tuple[bool, List[BlastHit]]:
    """
    Mark each hit as target or non-target.
    Returns (is_specific, annotated_hits).
    specific = True if ALL top hits contain the target organism name tokens.
    """
    target_tokens = set(target_organism.lower().split())

    for hit in hits:
        hit_tokens = set(hit.organism.lower().split())
        hit.is_target = bool(target_tokens & hit_tokens)

    non_target = [h for h in hits if not h.is_target]
    is_specific = len(non_target) == 0
    return is_specific, hits


# ── Synchronous single-sequence BLAST ─────────────────────────────────────────

def blast_sequence(
    sequence: str,
    label: str,
    target_organism: str,
    session: requests.Session = None,
) -> BlastResult:
    result = BlastResult(query_seq=sequence, query_label=label)
    sess   = session or requests.Session()

    logger.info("Submitting BLAST for %s (%s)", label, sequence)
    rid = _submit_blast(sequence, sess)
    if not rid:
        result.error = "Submit failed"
        return result

    logger.info("BLAST RID=%s — polling...", rid)
    xml = _poll_blast(rid, sess)
    if not xml:
        result.error = "Poll timeout or failure"
        return result

    hits = _parse_blast_xml(xml)
    specific, hits = _evaluate_specificity(hits, target_organism)
    result.hits     = hits
    result.specific = specific

    logger.info(
        "BLAST %s: %d hits, specific=%s",
        label, len(hits), specific
    )
    return result


# ── Parallel BLAST for a full PrimerSet ───────────────────────────────────────

def validate_primer_set(
    ps: PrimerSet,
    target_organism: str,
    workers: int = None,
) -> PrimerSet:
    """
    BLAST all three sequences in a PrimerSet in parallel using a thread pool.
    Updates ps.fwd_blast, ps.rev_blast, ps.probe_blast, ps.blast_pass in-place.
    """
    workers = workers or config.BLAST_WORKERS
    seqs = [
        (ps.fwd_seq,   "forward_primer"),
        (ps.rev_seq,   "reverse_primer"),
        (ps.probe_seq, "probe"),
    ]

    from concurrent.futures import ThreadPoolExecutor, as_completed
    results: Dict[str, BlastResult] = {}

    def _run(seq_label):
        seq, label = seq_label
        if not seq:
            return label, None
        sess = requests.Session()
        return label, blast_sequence(seq, label, target_organism, sess)

    with ThreadPoolExecutor(max_workers=min(workers, len(seqs))) as pool:
        futures = {pool.submit(_run, sl): sl[1] for sl in seqs}
        for fut in as_completed(futures):
            label, res = fut.result()
            results[label] = res

    ps.fwd_blast   = results.get("forward_primer")
    ps.rev_blast   = results.get("reverse_primer")
    ps.probe_blast = results.get("probe")

    # Overall pass: all three must be specific (or not have errors)
    checks = [
        r for r in [ps.fwd_blast, ps.rev_blast, ps.probe_blast]
        if r is not None and r.error is None
    ]
    ps.blast_pass = all(r.specific for r in checks) if checks else None

    return ps


def validate_all_sets(
    primer_sets: List[PrimerSet],
    target_organism: str,
    max_sets: int = 3,
) -> List[PrimerSet]:
    """
    BLAST the top *max_sets* primer sets (prioritise constraint-passing ones).
    Returns the list with blast fields populated.
    """
    to_blast = primer_sets[:max_sets]
    logger.info("BLASTing %d primer set(s) against '%s'", len(to_blast), target_organism)

    for i, ps in enumerate(to_blast):
        logger.info("--- BLAST run %d/%d ---", i + 1, len(to_blast))
        validate_primer_set(ps, target_organism)

    return primer_sets
