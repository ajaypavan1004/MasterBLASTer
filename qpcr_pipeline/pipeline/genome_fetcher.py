"""
genome_fetcher.py
Fetches reference genome sequences from NCBI for:
  - The target organism  (used as ROI source)
  - Close relatives      (used for specificity exclusion)
"""

import logging
import time
from typing import List, Optional

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from . import config

logger = logging.getLogger(__name__)


def _configure_entrez(email: str, api_key: str = ""):
    Entrez.email   = email or config.ENTREZ_EMAIL
    Entrez.tool    = config.ENTREZ_TOOL
    if api_key:
        Entrez.api_key = api_key


def search_taxon_accessions(
    organism: str,
    db: str = "nucleotide",
    max_seqs: int = 5,
    email: str = "",
    api_key: str = "",
) -> List[str]:
    """
    Search NCBI nucleotide DB for RefSeq complete genomes / sequences
    for *organism*.  Returns a list of accession IDs.
    """
    _configure_entrez(email, api_key)

    # Prefer RefSeq complete genomes; fall back to any sequence
    queries = [
        f'"{organism}"[Organism] AND "complete genome"[Title] AND refseq[filter]',
        f'"{organism}"[Organism] AND "complete sequence"[Title]',
        f'"{organism}"[Organism]',
    ]

    for q in queries:
        logger.info("Searching NCBI nucleotide: %s", q)
        try:
            handle = Entrez.esearch(db=db, term=q, retmax=max_seqs, usehistory="y")
            record = Entrez.read(handle)
            handle.close()
            ids = record.get("IdList", [])
            if ids:
                logger.info("Found %d accession(s) for '%s'", len(ids), organism)
                return ids
        except Exception as exc:
            logger.warning("NCBI search failed (%s): %s", q, exc)
            time.sleep(1)

    logger.warning("No NCBI records found for '%s'", organism)
    return []


def fetch_sequences(
    accession_ids: List[str],
    db: str = "nucleotide",
    email: str = "",
    api_key: str = "",
    batch_size: int = 5,
) -> List[SeqRecord]:
    """
    Download sequences by accession ID (GenBank format).
    Returns a list of SeqRecord objects.
    """
    _configure_entrez(email, api_key)
    records: List[SeqRecord] = []

    for i in range(0, len(accession_ids), batch_size):
        batch = accession_ids[i : i + batch_size]
        id_str = ",".join(batch)
        logger.info("Fetching sequences: %s", id_str)
        try:
            handle = Entrez.efetch(db=db, id=id_str, rettype="fasta", retmode="text")
            for rec in SeqIO.parse(handle, "fasta"):
                if len(rec.seq) > 0:
                    records.append(rec)
                    logger.info("  Fetched: %s  len=%d", rec.id, len(rec.seq))
            handle.close()
            time.sleep(0.4)          # be polite to NCBI
        except Exception as exc:
            logger.error("Fetch failed for IDs %s: %s", id_str, exc)
            time.sleep(2)

    return records


def get_close_relatives(
    organism: str,
    n: int = 3,
    email: str = "",
    api_key: str = "",
) -> List[SeqRecord]:
    """
    Heuristic: drop the species epithet and search at genus level,
    excluding the target species.  Used to build the exclusion set
    for ROI screening.
    """
    parts = organism.strip().split()
    if len(parts) < 2:
        return []

    genus = parts[0]
    species_query = " ".join(parts[:2])

    query = (
        f'"{genus}"[Organism] NOT "{species_query}"[Organism] '
        f'AND "complete genome"[Title]'
    )
    logger.info("Searching close relatives: %s", query)
    _configure_entrez(email, api_key)

    try:
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=n)
        record = Entrez.read(handle)
        handle.close()
        ids = record.get("IdList", [])
    except Exception as exc:
        logger.warning("Relative search failed: %s", exc)
        return []

    if not ids:
        return []

    return fetch_sequences(ids, email=email, api_key=api_key)
