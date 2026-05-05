"""
roi_selector.py
Identifies Regions of Interest (ROI) unique to the target organism
using k-mer based subtraction against close-relative sequences.

Strategy:
1. Build a k-mer set from all close-relative / exclusion sequences.
2. Slide a window across the target sequence; score each window by
   the fraction of k-mers NOT found in the exclusion set.
3. Rank windows by uniqueness score; return the top candidates.
"""

import logging
from dataclasses import dataclass
from typing import List, Optional, Tuple

from Bio.SeqRecord import SeqRecord

from . import config

logger = logging.getLogger(__name__)


@dataclass
class ROICandidate:
    record_id: str          # source sequence accession
    start: int              # 0-based start on source sequence
    end: int                # 0-based end (exclusive)
    sequence: str           # the actual nucleotide string
    uniqueness_score: float # fraction of k-mers absent from exclusion set (0–1)
    length: int

    def __str__(self):
        return (
            f"ROI [{self.record_id}:{self.start}-{self.end}] "
            f"len={self.length} uniqueness={self.uniqueness_score:.3f}"
        )


def _build_kmer_set(records: List[SeqRecord], k: int) -> set:
    """Build a set of all k-mers (and their reverse complements) from records."""
    kmer_set: set = set()
    for rec in records:
        seq = str(rec.seq).upper().replace("N", "")
        rc  = str(rec.seq.reverse_complement()).upper().replace("N", "")
        for s in (seq, rc):
            for i in range(len(s) - k + 1):
                kmer_set.add(s[i : i + k])
    logger.info("Built exclusion k-mer set: %d k-mers (k=%d)", len(kmer_set), k)
    return kmer_set


def _score_window(seq: str, exclusion_kmers: set, k: int) -> float:
    """
    Returns the fraction of k-mers in *seq* that are ABSENT from exclusion_kmers.
    Score of 1.0 = fully unique; 0.0 = fully shared with relatives.
    """
    total  = 0
    unique = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        total += 1
        if kmer not in exclusion_kmers:
            unique += 1
    return unique / total if total > 0 else 0.0


def find_roi_candidates(
    target_records: List[SeqRecord],
    exclusion_records: List[SeqRecord],
    window_size: int = 500,
    step: int = 100,
    k: int = None,
    min_uniqueness: float = 0.80,
    top_n: int = 10,
    min_length: int = None,
) -> List[ROICandidate]:
    """
    Slide a window over target sequences and score each window for uniqueness
    relative to *exclusion_records* (close relatives).

    Returns top_n ROICandidate objects sorted by uniqueness_score descending.
    """
    k          = k          or config.KMER_K
    min_length = min_length or config.ROI_MIN_LENGTH
    window_size = max(window_size, min_length)

    if not exclusion_records:
        logger.warning(
            "No exclusion sequences provided — ROI uniqueness cannot be verified. "
            "All windows will score 1.0 (assumed unique)."
        )
        exclusion_kmers: set = set()
    else:
        exclusion_kmers = _build_kmer_set(exclusion_records, k)

    candidates: List[ROICandidate] = []

    for rec in target_records:
        seq = str(rec.seq).upper()
        seqlen = len(seq)
        logger.info("Scanning %s (len=%d) with window=%d step=%d", rec.id, seqlen, window_size, step)

        for start in range(0, seqlen - window_size + 1, step):
            end    = start + window_size
            window = seq[start:end]

            # Skip windows with too many Ns
            n_frac = window.count("N") / len(window)
            if n_frac > 0.05:
                continue

            score = _score_window(window, exclusion_kmers, k)

            if score >= min_uniqueness:
                candidates.append(
                    ROICandidate(
                        record_id=rec.id,
                        start=start,
                        end=end,
                        sequence=window,
                        uniqueness_score=score,
                        length=window_size,
                    )
                )

    # Sort by score descending, then deduplicate overlapping windows
    candidates.sort(key=lambda c: c.uniqueness_score, reverse=True)
    candidates = _deduplicate(candidates, min_gap=window_size // 2)

    logger.info("Found %d ROI candidates (min_uniqueness=%.2f)", len(candidates), min_uniqueness)

    if not candidates:
        logger.warning(
            "No ROI candidates met uniqueness threshold %.2f. "
            "Try lowering --min-uniqueness.", min_uniqueness
        )

    return candidates[:top_n]


def _deduplicate(
    candidates: List[ROICandidate], min_gap: int = 250
) -> List[ROICandidate]:
    """
    Remove candidates that overlap heavily with a higher-scoring candidate
    from the same source record.
    """
    kept: List[ROICandidate] = []
    for cand in candidates:
        overlap = False
        for k_cand in kept:
            if k_cand.record_id != cand.record_id:
                continue
            # Check overlap
            overlap_start = max(cand.start, k_cand.start)
            overlap_end   = min(cand.end,   k_cand.end)
            if overlap_end - overlap_start > min_gap:
                overlap = True
                break
        if not overlap:
            kept.append(cand)
    return kept
