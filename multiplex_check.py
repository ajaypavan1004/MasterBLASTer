#!/usr/bin/env python3
"""
multiplex_check.py
==================
Cross-dimer and multiplex compatibility checker for qPCR panels.

Automatically finds the best combination of primer sets across all organisms
by trying all combinations of top-N sets per organism and scoring by:
  1. Fewest cross-dimer failures
  2. Smallest Tm spread
  3. Lowest total Primer3 penalty

Separate thresholds for primers (-6.0) and probes (-9.0) since probes
run at lower concentrations and probe-probe interactions are less critical.

Usage:
    python multiplex_check.py \\
        --panel "FAM:output/mycolicibacterium_fortuitum_results.json" \\
                "SUN:output/mycobacteroides_chelonae_results.json" \\
        --output output/multiplex_report.txt \\
        --top-n 3
"""

import argparse
import json
import logging
import os
import sys
from dataclasses import dataclass
from itertools import combinations, product
from typing import List, Tuple

import primer3

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("multiplex")

PRIMER_DIMER_THRESHOLD = -6.0   # kcal/mol — stricter for primers
PROBE_DIMER_THRESHOLD  = -9.0   # kcal/mol — more lenient for probes (lower conc.)
TM_SPREAD_MAX          = 5.0    # °C


# ── Data classes ──────────────────────────────────────────────────────────────

@dataclass
class PrimerEntry:
    channel:    str
    organism:   str
    rank:       int
    fwd_seq:    str
    rev_seq:    str
    probe_seq:  str
    fwd_tm:     float
    rev_tm:     float
    probe_tm:   float
    amplicon:   int
    penalty:    float


@dataclass
class DimerResult:
    label_a:   str
    label_b:   str
    dg:        float
    passes:    bool
    is_probe_pair: bool   # True if both sequences are probes


# ── Thermodynamics ────────────────────────────────────────────────────────────

def heterodimer_dg(seq_a: str, seq_b: str) -> float:
    try:
        return primer3.calc_heterodimer(seq_a, seq_b).dg / 1000.0
    except Exception:
        return 0.0

def homodimer_dg(seq: str) -> float:
    try:
        return primer3.calc_homodimer(seq).dg / 1000.0
    except Exception:
        return 0.0

def get_threshold(label_a: str, label_b: str) -> float:
    """Return appropriate threshold based on whether pair involves probes."""
    a_is_probe = label_a.endswith("_PRB")
    b_is_probe = label_b.endswith("_PRB")
    if a_is_probe and b_is_probe:
        return PROBE_DIMER_THRESHOLD   # probe-probe: lenient
    elif a_is_probe or b_is_probe:
        return (PRIMER_DIMER_THRESHOLD + PROBE_DIMER_THRESHOLD) / 2  # primer-probe: midpoint
    else:
        return PRIMER_DIMER_THRESHOLD  # primer-primer: strict


# ── Load candidates ───────────────────────────────────────────────────────────

def load_candidates(panel_args: List[str], top_n: int) -> List[List[PrimerEntry]]:
    all_candidates = []

    for arg in panel_args:
        if ":" not in arg:
            log.error("Panel entry must be 'CHANNEL:path' — got: %s", arg)
            sys.exit(1)
        channel, path = arg.split(":", 1)
        channel = channel.strip().upper()

        if not os.path.exists(path):
            log.error("File not found: %s", path)
            sys.exit(1)

        with open(path) as f:
            data = json.load(f)

        organism = data.get("target_organism", "Unknown")
        sets = data.get("primer_sets", [])

        usable = [ps for ps in sets if ps.get("passed_constraints") and ps.get("blast", {}).get("overall_pass")]
        if not usable:
            usable = [ps for ps in sets if ps.get("passed_constraints")]
        if not usable:
            usable = sets

        candidates = []
        for rank, ps in enumerate(usable[:top_n], start=1):
            primers = ps.get("primers", {})
            candidates.append(PrimerEntry(
                channel   = channel,
                organism  = organism,
                rank      = rank,
                fwd_seq   = primers["forward"]["sequence"],
                rev_seq   = primers["reverse"]["sequence"],
                probe_seq = primers["probe"]["sequence"],
                fwd_tm    = primers["forward"]["tm"],
                rev_tm    = primers["reverse"]["tm"],
                probe_tm  = primers["probe"]["tm"],
                amplicon  = ps.get("amplicon_size", 0),
                penalty   = ps.get("primer3_penalty", 0),
            ))

        log.info("Loaded %s (%s) — %d candidate sets", channel, organism, len(candidates))
        all_candidates.append(candidates)

    return all_candidates


# ── Scoring ───────────────────────────────────────────────────────────────────

def score_combination(combo: List[PrimerEntry]) -> Tuple[int, float, float]:
    seqs = []
    for e in combo:
        seqs.append((f"{e.channel}_FWD", e.fwd_seq))
        seqs.append((f"{e.channel}_REV", e.rev_seq))
        if e.probe_seq:
            seqs.append((f"{e.channel}_PRB", e.probe_seq))

    failures = 0
    for (la, sa), (lb, sb) in combinations(seqs, 2):
        thresh = get_threshold(la, lb)
        if heterodimer_dg(sa, sb) < thresh:
            failures += 1

    tms = [e.fwd_tm for e in combo] + [e.rev_tm for e in combo]
    spread = max(tms) - min(tms)
    total_penalty = sum(e.penalty for e in combo)

    return failures, spread, total_penalty


def find_best_combination(all_candidates: List[List[PrimerEntry]]) -> Tuple[List[PrimerEntry], int, float]:
    total_combos = 1
    for candidates in all_candidates:
        total_combos *= len(candidates)
    log.info("Trying %d combinations...", total_combos)

    best_combo = None
    best_score = (999, 999.0, 999.0)

    for combo in product(*all_candidates):
        combo = list(combo)
        score = score_combination(combo)
        if score < best_score:
            best_score = score
            best_combo = combo

    return best_combo, best_score[0], best_score[1]


# ── Cross-dimer matrix ────────────────────────────────────────────────────────

def run_cross_dimer_matrix(entries: List[PrimerEntry]) -> List[DimerResult]:
    seqs = []
    for e in entries:
        seqs.append((f"{e.channel}_FWD", e.fwd_seq))
        seqs.append((f"{e.channel}_REV", e.rev_seq))
        if e.probe_seq:
            seqs.append((f"{e.channel}_PRB", e.probe_seq))

    results = []
    for (la, sa), (lb, sb) in combinations(seqs, 2):
        dg = heterodimer_dg(sa, sb)
        thresh = get_threshold(la, lb)
        is_probe_pair = la.endswith("_PRB") and lb.endswith("_PRB")
        results.append(DimerResult(
            label_a=la, label_b=lb, dg=dg,
            passes=dg >= thresh,
            is_probe_pair=is_probe_pair,
        ))
    return results


# ── Tm spread ─────────────────────────────────────────────────────────────────

def check_tm_spread(entries: List[PrimerEntry]) -> Tuple[bool, float, float, float]:
    tms = [e.fwd_tm for e in entries] + [e.rev_tm for e in entries]
    min_tm, max_tm = min(tms), max(tms)
    spread = max_tm - min_tm
    return spread <= TM_SPREAD_MAX, min_tm, max_tm, spread


# ── Report ────────────────────────────────────────────────────────────────────

def write_report(entries: List[PrimerEntry], dimer_results: List[DimerResult],
                 output_path: str, searched_combos: int):
    lines = []
    lines.append("=" * 70)
    lines.append("  qPCR Multiplex Compatibility Report")
    lines.append("=" * 70)
    lines.append(f"  Searched {searched_combos} combinations — showing best result")
    lines.append(f"  Thresholds: primers ≥ {PRIMER_DIMER_THRESHOLD} kcal/mol | "
                 f"probes ≥ {PROBE_DIMER_THRESHOLD} kcal/mol "
                 f"(probes run at lower concentration — less critical)")
    lines.append("")

    lines.append("BEST COMBINATION — PANEL SUMMARY")
    lines.append("─" * 70)
    for e in entries:
        mean_tm = (e.fwd_tm + e.rev_tm) / 2
        lines.append(f"  {e.channel:<8} {e.organism}  (Rank #{e.rank})")
        lines.append(f"           FWD  5'-{e.fwd_seq}-3'  Tm={e.fwd_tm:.1f}°C  ΔG={homodimer_dg(e.fwd_seq):.1f}")
        lines.append(f"           REV  5'-{e.rev_seq}-3'  Tm={e.rev_tm:.1f}°C  ΔG={homodimer_dg(e.rev_seq):.1f}")
        if e.probe_seq:
            lines.append(f"           PRB  5'-{e.probe_seq}-3'  Tm={e.probe_tm:.1f}°C  ΔG={homodimer_dg(e.probe_seq):.1f}")
        lines.append(f"           Amplicon={e.amplicon}nt  Mean primer Tm={mean_tm:.1f}°C")
        lines.append("")

    lines.append("TM SPREAD CHECK")
    lines.append("─" * 70)
    tm_ok, min_tm, max_tm, spread = check_tm_spread(entries)
    lines.append(f"  Primer Tm range : {min_tm:.1f}°C – {max_tm:.1f}°C")
    lines.append(f"  Spread          : {spread:.1f}°C  (threshold ≤{TM_SPREAD_MAX}°C)  {'✓ PASS' if tm_ok else '✗ FAIL'}")
    lines.append("")
    for e in entries:
        lines.append(f"  {e.channel:<8}  FWD={e.fwd_tm:.1f}°C  REV={e.rev_tm:.1f}°C  PRB={e.probe_tm:.1f}°C")
    lines.append("")

    lines.append("CROSS-DIMER MATRIX")
    lines.append("─" * 70)
    failures = [r for r in dimer_results if not r.passes]
    probe_failures = [r for r in failures if r.is_probe_pair]
    primer_failures = [r for r in failures if not r.is_probe_pair]
    lines.append(f"  Total pairs checked : {len(dimer_results)}")
    lines.append(f"  Passed              : {len(dimer_results) - len(failures)}")
    lines.append(f"  Failed              : {len(failures)}  "
                 f"(primer-primer/primer-probe: {len(primer_failures)} | probe-probe: {len(probe_failures)})")
    lines.append("")

    if failures:
        lines.append("  ✗ PROBLEMATIC PAIRS:")
        for r in sorted(failures, key=lambda x: x.dg):
            tag = "[probe-probe]" if r.is_probe_pair else "[primer]"
            thresh = get_threshold(r.label_a, r.label_b)
            lines.append(f"    {r.label_a:<15} × {r.label_b:<15}  ΔG={r.dg:.2f}  threshold={thresh:.1f}  {tag}")
        lines.append("")
    else:
        lines.append("  ✓ No problematic cross-dimers detected.")
        lines.append("")

    lines.append("  Full matrix:")
    lines.append(f"  {'Pair A':<15} × {'Pair B':<15}  ΔG (kcal/mol)  Threshold  Status")
    lines.append(f"  {'─'*15}   {'─'*15}  {'─'*13}  {'─'*9}  {'─'*6}")
    for r in sorted(dimer_results, key=lambda x: x.dg):
        thresh = get_threshold(r.label_a, r.label_b)
        lines.append(f"  {r.label_a:<15} × {r.label_b:<15}  {r.dg:>8.2f}       {thresh:>6.1f}     {'✓' if r.passes else '✗'}")
    lines.append("")

    dimer_ok = len(failures) == 0
    overall  = tm_ok and dimer_ok
    lines.append("MULTIPLEX COMPATIBILITY VERDICT")
    lines.append("─" * 70)
    lines.append(f"  Tm spread         : {'✓ PASS' if tm_ok else '✗ FAIL'}")
    lines.append(f"  Primer cross-dimers: {'✓ PASS' if len(primer_failures) == 0 else f'✗ FAIL ({len(primer_failures)} pairs)'}")
    lines.append(f"  Probe cross-dimers : {'✓ PASS' if len(probe_failures) == 0 else f'✗ FAIL ({len(probe_failures)} pairs)'}")
    lines.append(f"  Overall            : {'✓ COMPATIBLE — safe to multiplex' if overall else '✗ NOT COMPATIBLE — review failures above'}")
    lines.append("")

    lines.append("CHANNEL ASSIGNMENT")
    lines.append("─" * 70)
    lines.append(f"  {'Channel':<10} {'Organism':<40} {'Amplicon'}")
    lines.append(f"  {'─'*10} {'─'*40} {'─'*8}")
    for e in entries:
        lines.append(f"  {e.channel:<10} {e.organism:<40} {e.amplicon}nt")
    lines.append("")

    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else ".", exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write("\n".join(lines))
    log.info("Multiplex report → %s", output_path)
    print("\n" + "\n".join(lines))


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    global PRIMER_DIMER_THRESHOLD, PROBE_DIMER_THRESHOLD, TM_SPREAD_MAX
    p = argparse.ArgumentParser(description="Multiplex qPCR compatibility checker")
    p.add_argument("--panel", nargs="+", required=True, metavar="CHANNEL:JSON_PATH")
    p.add_argument("--output", default="output/multiplex_report.txt")
    p.add_argument("--top-n", type=int, default=3,
                   help="Number of top sets per organism to try (default: 3)")
    p.add_argument("--dimer-threshold", type=float, default=PRIMER_DIMER_THRESHOLD,
                   help=f"Primer heterodimer ΔG threshold (default: {PRIMER_DIMER_THRESHOLD})")
    p.add_argument("--probe-dimer-threshold", type=float, default=PROBE_DIMER_THRESHOLD,
                   help=f"Probe heterodimer ΔG threshold (default: {PROBE_DIMER_THRESHOLD})")
    p.add_argument("--tm-spread", type=float, default=TM_SPREAD_MAX,
                   help=f"Max Tm spread in °C (default: {TM_SPREAD_MAX})")
    args = p.parse_args()

    PRIMER_DIMER_THRESHOLD = args.dimer_threshold
    PROBE_DIMER_THRESHOLD  = args.probe_dimer_threshold
    TM_SPREAD_MAX          = args.tm_spread

    log.info("Loading panel (top-%d sets per organism)...", args.top_n)
    all_candidates = load_candidates(args.panel, args.top_n)

    total_combos = 1
    for c in all_candidates: total_combos *= len(c)

    best_combo, n_failures, tm_spread = find_best_combination(all_candidates)

    log.info("Best combination: %d failures, %.1f°C Tm spread", n_failures, tm_spread)
    for e in best_combo:
        log.info("  %s → Rank #%d  %s", e.channel, e.rank, e.organism)

    dimer_results = run_cross_dimer_matrix(best_combo)
    write_report(best_combo, dimer_results, args.output, total_combos)


if __name__ == "__main__":
    main()