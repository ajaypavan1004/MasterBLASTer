"""
reporter.py
Generates pipeline output:
  - JSON  (machine-readable, full detail)
  - TSV   (spreadsheet-friendly summary)
  - TXT   (human-readable report)
"""

import json
import logging
import os
from datetime import datetime
from typing import List, Optional

from .primer_designer import PrimerSet
from . import config

logger = logging.getLogger(__name__)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _blast_summary(blast_result) -> dict:
    if blast_result is None:
        return {"status": "not_run"}
    if blast_result.error:
        return {"status": "error", "message": blast_result.error}
    return {
        "status":   "ok",
        "specific": blast_result.specific,
        "top_hits": [
            {
                "rank":       h.rank,
                "accession":  h.accession,
                "organism":   h.organism,
                "identity":   round(h.identity, 1),
                "coverage":   round(h.coverage, 1),
                "evalue":     h.evalue,
                "is_target":  h.is_target,
            }
            for h in blast_result.hits
        ],
    }


def _set_to_dict(ps: PrimerSet, target_organism: str) -> dict:
    return {
        "pair_index":           ps.pair_index,
        "roi": {
            "source_id":        ps.roi.record_id,
            "start":            ps.roi.start,
            "end":              ps.roi.end,
            "uniqueness_score": round(ps.roi.uniqueness_score, 4),
        },
        "primers": {
            "forward": {
                "sequence": ps.fwd_seq,
                "tm":       round(ps.fwd_tm, 2),
                "gc":       round(ps.fwd_gc, 1),
                "length":   len(ps.fwd_seq),
                "start":    ps.fwd_start,
            },
            "reverse": {
                "sequence": ps.rev_seq,
                "tm":       round(ps.rev_tm, 2),
                "gc":       round(ps.rev_gc, 1),
                "length":   len(ps.rev_seq),
                "start":    ps.rev_start,
            },
            "probe": {
                "sequence": ps.probe_seq,
                "tm":       round(ps.probe_tm, 2),
                "gc":       round(ps.probe_gc, 1),
                "length":   len(ps.probe_seq),
                "start":    ps.probe_start,
            },
        },
        "amplicon_size":        ps.amplicon_size,
        "primer3_penalty":      round(ps.primer3_penalty, 4),
        "passed_constraints":   ps.passed_constraints,
        "constraint_notes":     ps.constraint_notes,
        "blast": {
            "forward_primer":   _blast_summary(ps.fwd_blast),
            "reverse_primer":   _blast_summary(ps.rev_blast),
            "probe":            _blast_summary(ps.probe_blast),
            "overall_pass":     ps.blast_pass,
        },
    }


# ── Public API ────────────────────────────────────────────────────────────────

def write_json(
    primer_sets: List[PrimerSet],
    target_organism: str,
    output_dir: str = None,
    filename: str = None,
) -> str:
    output_dir = output_dir or config.OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    filename = filename or f"{_slug(target_organism)}_results.json"
    path = os.path.join(output_dir, filename)

    payload = {
        "pipeline":       "qpcr_primer_pipeline",
        "version":        "1.0.0",
        "timestamp":      datetime.utcnow().isoformat() + "Z",
        "target_organism": target_organism,
        "n_sets":         len(primer_sets),
        "primer_sets":    [_set_to_dict(ps, target_organism) for ps in primer_sets],
    }

    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)

    logger.info("JSON report → %s", path)
    return path


def write_tsv(
    primer_sets: List[PrimerSet],
    target_organism: str,
    output_dir: str = None,
    filename: str = None,
) -> str:
    output_dir = output_dir or config.OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    filename = filename or f"{_slug(target_organism)}_results.tsv"
    path = os.path.join(output_dir, filename)

    headers = [
        "pair_index", "roi_id", "roi_start", "roi_end", "roi_uniqueness",
        "fwd_seq", "fwd_len", "fwd_tm", "fwd_gc",
        "rev_seq", "rev_len", "rev_tm", "rev_gc",
        "probe_seq", "probe_len", "probe_tm", "probe_gc",
        "amplicon_size", "primer3_penalty",
        "passed_constraints", "constraint_notes",
        "blast_fwd_specific", "blast_rev_specific", "blast_probe_specific",
        "blast_overall_pass",
    ]

    rows = []
    for ps in primer_sets:
        def spec(br):
            if br is None: return "N/A"
            if br.error:   return f"ERROR:{br.error}"
            return str(br.specific)

        rows.append([
            ps.pair_index,
            ps.roi.record_id, ps.roi.start, ps.roi.end,
            round(ps.roi.uniqueness_score, 4),
            ps.fwd_seq, len(ps.fwd_seq), round(ps.fwd_tm, 2), round(ps.fwd_gc, 1),
            ps.rev_seq, len(ps.rev_seq), round(ps.rev_tm, 2), round(ps.rev_gc, 1),
            ps.probe_seq, len(ps.probe_seq), round(ps.probe_tm, 2), round(ps.probe_gc, 1),
            ps.amplicon_size, round(ps.primer3_penalty, 4),
            ps.passed_constraints, "; ".join(ps.constraint_notes),
            spec(ps.fwd_blast), spec(ps.rev_blast), spec(ps.probe_blast),
            str(ps.blast_pass),
        ])

    with open(path, "w") as fh:
        fh.write("\t".join(headers) + "\n")
        for row in rows:
            fh.write("\t".join(str(v) for v in row) + "\n")

    logger.info("TSV report → %s", path)
    return path


def write_txt(
    primer_sets: List[PrimerSet],
    target_organism: str,
    output_dir: str = None,
    filename: str = None,
) -> str:
    output_dir = output_dir or config.OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    filename = filename or f"{_slug(target_organism)}_report.txt"
    path = os.path.join(output_dir, filename)

    lines = [
        "=" * 70,
        f"  qPCR Primer/Probe Pipeline Report",
        f"  Target: {target_organism}",
        f"  Generated: {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 70,
        "",
    ]

    passed = [ps for ps in primer_sets if ps.passed_constraints and ps.blast_pass]
    lines.append(f"Total primer sets designed : {len(primer_sets)}")
    lines.append(f"Passed ALL constraints     : {len(passed)}")
    lines.append(f"BLAST-validated            : {sum(1 for ps in primer_sets if ps.blast_pass)}")
    lines.append("")

    for ps in primer_sets:
        blast_ok = ps.blast_pass if ps.blast_pass is not None else True
        status = "✓ PASS" if (ps.passed_constraints and blast_ok) else "✗ FAIL"
        lines.append(f"{'─'*70}")
        lines.append(f"  Pair #{ps.pair_index}  {status}  |  ROI {ps.roi.record_id}:{ps.roi.start}-{ps.roi.end}")
        lines.append(f"  ROI uniqueness score: {ps.roi.uniqueness_score:.3f}")
        lines.append("")
        lines.append(f"  FWD  5'-{ps.fwd_seq}-3'")
        lines.append(f"       Tm={ps.fwd_tm:.1f}°C  GC={ps.fwd_gc:.1f}%  len={len(ps.fwd_seq)}nt")
        lines.append("")
        lines.append(f"  REV  5'-{ps.rev_seq}-3'")
        lines.append(f"       Tm={ps.rev_tm:.1f}°C  GC={ps.rev_gc:.1f}%  len={len(ps.rev_seq)}nt")
        lines.append("")
        lines.append(f"  PRB  5'-{ps.probe_seq}-3'")
        lines.append(f"       Tm={ps.probe_tm:.1f}°C  GC={ps.probe_gc:.1f}%  len={len(ps.probe_seq)}nt")
        lines.append(f"       ΔTm vs primers = {ps.probe_tm - (ps.fwd_tm+ps.rev_tm)/2:+.1f}°C")
        lines.append("")
        lines.append(f"  Amplicon size : {ps.amplicon_size} nt")
        lines.append(f"  P3 penalty    : {ps.primer3_penalty:.4f}")

        if ps.constraint_notes:
            lines.append(f"  Constraint issues:")
            for note in ps.constraint_notes:
                lines.append(f"    • {note}")

        if ps.blast_pass is not None:
            lines.append(f"  BLAST overall : {'SPECIFIC ✓' if ps.blast_pass else 'NOT SPECIFIC ✗'}")
            for label, br in [("FWD", ps.fwd_blast), ("REV", ps.rev_blast), ("PRB", ps.probe_blast)]:
                if br and not br.error and br.hits:
                    top3 = br.hits[:3]
                    lines.append(f"    {label} top hits:")
                    for h in top3:
                        flag = "✓" if h.is_target else "✗"
                        lines.append(
                            f"      {flag} [{h.rank}] {h.organism:<30} "
                            f"id={h.identity:.1f}% e={h.evalue:.1e}"
                        )
        lines.append("")

    with open(path, "w") as fh:
        fh.write("\n".join(lines))

    logger.info("TXT report → %s", path)
    return path


def _slug(s: str) -> str:
    return s.lower().replace(" ", "_").replace("/", "_")[:40]
