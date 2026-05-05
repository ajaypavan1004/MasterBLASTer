"""
primer_designer.py
Wraps primer3-py to design forward primer, reverse primer, and hydrolysis probe
from an ROI sequence.  All constraints are pulled from config.py and can be
overridden per-call.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any

import primer3

from .roi_selector import ROICandidate
from . import config

logger = logging.getLogger(__name__)


@dataclass
class PrimerSet:
    roi: ROICandidate
    pair_index: int                     # Primer3 pair rank (0 = best)

    fwd_seq:    str = ""
    rev_seq:    str = ""
    probe_seq:  str = ""

    fwd_tm:     float = 0.0
    rev_tm:     float = 0.0
    probe_tm:   float = 0.0

    fwd_gc:     float = 0.0
    rev_gc:     float = 0.0
    probe_gc:   float = 0.0

    fwd_start:  int = 0                 # 0-based on ROI
    rev_start:  int = 0
    probe_start:int = 0

    amplicon_size: int = 0

    primer3_penalty: float = 0.0        # lower = better
    passed_constraints: bool = False
    constraint_notes: List[str] = field(default_factory=list)

    # BLAST results filled in later
    fwd_blast:   Optional[Dict] = None
    rev_blast:   Optional[Dict] = None
    probe_blast: Optional[Dict] = None
    blast_pass:  Optional[bool] = None

    def sequences(self) -> Dict[str, str]:
        return {
            "forward_primer": self.fwd_seq,
            "reverse_primer": self.rev_seq,
            "probe":          self.probe_seq,
        }

    def __str__(self):
        return (
            f"PrimerSet #{self.pair_index} | ROI {self.roi.record_id}:{self.roi.start}-{self.roi.end}\n"
            f"  FWD  {self.fwd_seq:<28} Tm={self.fwd_tm:.1f}°C GC={self.fwd_gc:.1f}%\n"
            f"  REV  {self.rev_seq:<28} Tm={self.rev_tm:.1f}°C GC={self.rev_gc:.1f}%\n"
            f"  PRB  {self.probe_seq:<28} Tm={self.probe_tm:.1f}°C GC={self.probe_gc:.1f}%\n"
            f"  Amplicon={self.amplicon_size} nt  Penalty={self.primer3_penalty:.2f}  "
            f"Pass={self.passed_constraints}"
        )


def _build_primer3_globals(overrides: Dict[str, Any] = None) -> Dict[str, Any]:
    """Build the PRIMER_GLOBALS dict for primer3."""
    cfg = {
        # Primer sizing
        "PRIMER_MIN_SIZE":       config.PRIMER_MIN_SIZE,
        "PRIMER_OPT_SIZE":       config.PRIMER_OPT_SIZE,
        "PRIMER_MAX_SIZE":       config.PRIMER_MAX_SIZE,
        # Primer Tm
        "PRIMER_MIN_TM":         config.PRIMER_MIN_TM,
        "PRIMER_OPT_TM":         config.PRIMER_OPT_TM,
        "PRIMER_MAX_TM":         config.PRIMER_MAX_TM,
        # Primer GC
        "PRIMER_MIN_GC":         config.PRIMER_MIN_GC,
        "PRIMER_OPT_GC":         config.PRIMER_OPT_GC,
        "PRIMER_MAX_GC":         config.PRIMER_MAX_GC,
        # Amplicon
        "PRIMER_PRODUCT_SIZE_RANGE": [[config.AMPLICON_MIN, config.AMPLICON_MAX]],
        # Internal oligo (probe)
        "PRIMER_INTERNAL_MIN_SIZE":  config.PROBE_MIN_SIZE,
        "PRIMER_INTERNAL_OPT_SIZE":  config.PROBE_OPT_SIZE,
        "PRIMER_INTERNAL_MAX_SIZE":  config.PROBE_MAX_SIZE,
        "PRIMER_INTERNAL_MIN_GC":    config.PROBE_MIN_GC,
        "PRIMER_INTERNAL_MAX_GC":    config.PROBE_MAX_GC,
        # Probe Tm — aim 7°C above primer opt Tm so the delta lands in [5,10]
        "PRIMER_INTERNAL_MIN_TM":    config.PRIMER_OPT_TM + config.PROBE_TM_DELTA_MIN,
        "PRIMER_INTERNAL_OPT_TM":    config.PRIMER_OPT_TM + (config.PROBE_TM_DELTA_MIN + config.PROBE_TM_DELTA_MAX) / 2,
        "PRIMER_INTERNAL_MAX_TM":    config.PRIMER_MAX_TM + config.PROBE_TM_DELTA_MAX,
        # General
        "PRIMER_NUM_RETURN":         10,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_MAX_SELF_ANY":       8,
        "PRIMER_MAX_SELF_END":       3,
        "PRIMER_PAIR_MAX_COMPL_ANY": 8,
        "PRIMER_PAIR_MAX_COMPL_END": 3,
        "PRIMER_MAX_POLY_X":         4,
        "PRIMER_INTERNAL_MAX_SELF_ANY": 8,
    }
    if overrides:
        cfg.update(overrides)
    return cfg


def _validate_set(ps: PrimerSet) -> PrimerSet:
    """
    Apply post-Primer3 constraint checks (especially probe Tm delta).
    Populates ps.passed_constraints and ps.constraint_notes.
    """
    notes = []
    ok    = True

    # Primer Tm range
    for name, tm in [("FWD", ps.fwd_tm), ("REV", ps.rev_tm)]:
        if not (config.PRIMER_MIN_TM <= tm <= config.PRIMER_MAX_TM):
            notes.append(f"{name} Tm {tm:.1f}°C outside [{config.PRIMER_MIN_TM},{config.PRIMER_MAX_TM}]")
            ok = False

    # Probe Tm delta above mean primer Tm
    mean_primer_tm = (ps.fwd_tm + ps.rev_tm) / 2
    delta = ps.probe_tm - mean_primer_tm
    if not (config.PROBE_TM_DELTA_MIN <= delta <= config.PROBE_TM_DELTA_MAX):
        notes.append(
            f"Probe Tm delta {delta:.1f}°C outside [{config.PROBE_TM_DELTA_MIN},{config.PROBE_TM_DELTA_MAX}]"
        )
        ok = False

    # Amplicon size
    if not (config.AMPLICON_MIN <= ps.amplicon_size <= config.AMPLICON_MAX):
        notes.append(f"Amplicon {ps.amplicon_size} nt outside [{config.AMPLICON_MIN},{config.AMPLICON_MAX}]")
        ok = False

    ps.passed_constraints = ok
    ps.constraint_notes   = notes
    return ps


def design_primers(
    roi: ROICandidate,
    n_pairs: int = 5,
    primer3_overrides: Dict[str, Any] = None,
) -> List[PrimerSet]:
    """
    Run Primer3 on an ROI and return up to *n_pairs* PrimerSet objects,
    sorted by primer3_penalty (ascending = better).

    Primer3 is asked to design BOTH flanking primers AND an internal
    probe (PRIMER_PICK_INTERNAL_OLIGO=1).
    """
    seq = roi.sequence.upper()
    seq_args = {
        "SEQUENCE_ID":       f"{roi.record_id}_{roi.start}_{roi.end}",
        "SEQUENCE_TEMPLATE": seq,
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [],  # no hard exclusion
    }

    global_args = _build_primer3_globals(primer3_overrides)
    global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 1   # design probe
    global_args["PRIMER_NUM_RETURN"]           = n_pairs

    logger.info("Running Primer3 on ROI %s:%d-%d (len=%d)",
                roi.record_id, roi.start, roi.end, roi.length)

    try:
        result = primer3.bindings.design_primers(seq_args, global_args)
    except Exception as exc:
        logger.error("Primer3 failed: %s", exc)
        return []

    n_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    logger.info("Primer3 returned %d pair(s)", n_returned)

    sets: List[PrimerSet] = []
    for i in range(n_returned):
        try:
            ps = PrimerSet(
                roi=roi,
                pair_index=i,
                fwd_seq=result[f"PRIMER_LEFT_{i}_SEQUENCE"],
                rev_seq=result[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                probe_seq=result.get(f"PRIMER_INTERNAL_{i}_SEQUENCE", ""),
                fwd_tm=result[f"PRIMER_LEFT_{i}_TM"],
                rev_tm=result[f"PRIMER_RIGHT_{i}_TM"],
                probe_tm=result.get(f"PRIMER_INTERNAL_{i}_TM", 0.0),
                fwd_gc=result[f"PRIMER_LEFT_{i}_GC_PERCENT"],
                rev_gc=result[f"PRIMER_RIGHT_{i}_GC_PERCENT"],
                probe_gc=result.get(f"PRIMER_INTERNAL_{i}_GC_PERCENT", 0.0),
                fwd_start=result[f"PRIMER_LEFT_{i}"][0],
                rev_start=result[f"PRIMER_RIGHT_{i}"][0],
                probe_start=result.get(f"PRIMER_INTERNAL_{i}", [0])[0],
                amplicon_size=result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                primer3_penalty=result[f"PRIMER_PAIR_{i}_PENALTY"],
            )
            ps = _validate_set(ps)
            sets.append(ps)
        except KeyError as exc:
            logger.debug("Missing Primer3 key for pair %d: %s", i, exc)

    # Sort: passing sets first, then by penalty
    sets.sort(key=lambda p: (not p.passed_constraints, p.primer3_penalty))
    return sets
