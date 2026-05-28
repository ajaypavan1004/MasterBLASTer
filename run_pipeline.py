#!/usr/bin/env python3
"""
qpcr_pipeline — CLI entry point
================================
Designs and validates qPCR primer/probe sets for any target organism.

Usage examples
--------------
# Minimal (uses NCBI API, email required):
python run_pipeline.py \
    --organism "Cyclospora cayetanensis" \
    --email your@email.com

# With API key (higher rate limits):
python run_pipeline.py \
    --organism "Gardnerella vaginalis" \
    --email your@email.com \
    --api-key YOUR_NCBI_API_KEY

# Relax primer length for difficult organisms:
python run_pipeline.py \
    --organism "Cyclospora cayetanensis" \
    --email your@email.com \
    --primer-max-size 26

# Skip BLAST (design-only mode):
python run_pipeline.py \
    --organism "Cyclospora cayetanensis" \
    --email your@email.com \
    --no-blast

# Use pre-downloaded FASTA instead of fetching from NCBI:
python run_pipeline.py \
    --organism "Cyclospora cayetanensis" \
    --email your@email.com \
    --target-fasta target.fasta \
    --exclusion-fasta relatives.fasta
"""

import argparse
import logging
import sys
import os
from typing import List

from Bio import SeqIO

# ── Pipeline modules ──────────────────────────────────────────────────────────
from pipeline import config
from pipeline.genome_fetcher   import search_taxon_accessions, fetch_sequences, get_close_relatives, search_gene_sequences, get_close_relative_gene_sequences
from pipeline.roi_selector     import find_roi_candidates, ROICandidate
from pipeline.primer_designer  import design_primers, PrimerSet
from pipeline.blast_validator  import validate_all_sets
from pipeline.reporter         import write_json, write_tsv, write_txt


# ── Logging setup ─────────────────────────────────────────────────────────────

def _setup_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s  %(levelname)-8s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )


# ── Argument parsing ──────────────────────────────────────────────────────────

def _parse_args():
    p = argparse.ArgumentParser(
        description="qPCR Primer/Probe Design Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Target
    p.add_argument("--organism", required=True,
                   help='Target organism (e.g. "Cyclospora cayetanensis")')
    p.add_argument("--email",    required=True,
                   help="Email address for NCBI Entrez (required by NCBI)")
    p.add_argument("--api-key",  default="",
                   help="NCBI API key (optional; increases rate limits)")

    # Input overrides
    p.add_argument("--target-fasta",    default="",
                   help="Path to pre-downloaded target FASTA (skips NCBI fetch)")
    p.add_argument("--exclusion-fasta", default="",
                   help="Path to exclusion/relatives FASTA (skips NCBI fetch)")
    p.add_argument("--max-target-seqs", type=int, default=3,
                   help="Max target sequences to fetch from NCBI (default: 3)")
    p.add_argument("--max-relative-seqs", type=int, default=10,
                   help="Max close-relative sequences to fetch (default: 10)")
    p.add_argument("--extra-exclusion", nargs="+", default=[],
                   metavar="ACCESSION",
                   help="Extra NCBI accessions to add to exclusion set (e.g. CP019058.1)")
    p.add_argument("--target-gene", default="",
                   help="Target a specific gene instead of whole genome (e.g. cpn60, rpoB, rnaseh). "
                        "Fetches gene sequences for ROI design and uses relative gene seqs for exclusion.")

    # ROI options
    p.add_argument("--roi-window",      type=int, default=500,
                   help="ROI window size in bp (default: 500)")
    p.add_argument("--roi-step",        type=int, default=100,
                   help="ROI sliding step in bp (default: 100)")
    p.add_argument("--min-uniqueness",  type=float, default=0.80,
                   help="Min k-mer uniqueness score 0-1 (default: 0.80)")
    p.add_argument("--top-rois",        type=int, default=5,
                   help="Number of top ROIs to pass to Primer3 (default: 5)")

    # Primer constraints
    p.add_argument("--primer-min-size", type=int,   default=config.PRIMER_MIN_SIZE)
    p.add_argument("--primer-opt-size", type=int,   default=config.PRIMER_OPT_SIZE)
    p.add_argument("--primer-max-size", type=int,   default=config.PRIMER_MAX_SIZE)
    p.add_argument("--primer-min-tm",   type=float, default=config.PRIMER_MIN_TM)
    p.add_argument("--primer-opt-tm",   type=float, default=config.PRIMER_OPT_TM)
    p.add_argument("--primer-max-tm",   type=float, default=config.PRIMER_MAX_TM)

    # BLAST
    p.add_argument("--no-blast", action="store_true",
                   help="Skip BLAST validation (design-only mode)")
    p.add_argument("--blast-sets", type=int, default=3,
                   help="Number of top primer sets to BLAST (default: 3)")
    p.add_argument("--blast-workers", type=int, default=config.BLAST_WORKERS,
                   help=f"Parallel BLAST threads (default: {config.BLAST_WORKERS})")

    # Output
    p.add_argument("--output-dir", default="output",
                   help="Directory for output files (default: ./output)")
    p.add_argument("--verbose", "-v", action="store_true")

    return p.parse_args()


# ── Pipeline stages ───────────────────────────────────────────────────────────

def stage_fetch(args) -> tuple:
    """Stage 1: Fetch or load target and exclusion sequences."""
    log = logging.getLogger("stage.fetch")

    if args.target_fasta:
        log.info("Loading target sequences from %s", args.target_fasta)
        target_recs = list(SeqIO.parse(args.target_fasta, "fasta"))
    elif args.target_gene:
        log.info("Gene-targeted mode: fetching '%s' sequences for '%s'", args.target_gene, args.organism)
        target_recs = search_gene_sequences(
            args.organism,
            args.target_gene,
            max_seqs=args.max_target_seqs,
            email=args.email,
            api_key=args.api_key,
        )
    else:
        log.info("Fetching target sequences from NCBI: %s", args.organism)
        ids = search_taxon_accessions(
            args.organism,
            max_seqs=args.max_target_seqs,
            email=args.email,
            api_key=args.api_key,
        )
        if not ids:
            log.error("No sequences found for '%s'. Exiting.", args.organism)
            sys.exit(1)
        target_recs = fetch_sequences(ids, email=args.email, api_key=args.api_key)

    if not target_recs:
        log.error("Failed to load any target sequences. Exiting.")
        sys.exit(1)

    log.info("Loaded %d target sequence(s)", len(target_recs))

    if args.exclusion_fasta:
        log.info("Loading exclusion sequences from %s", args.exclusion_fasta)
        excl_recs = list(SeqIO.parse(args.exclusion_fasta, "fasta"))
    elif args.target_gene:
        log.info("Gene-targeted mode: fetching '%s' sequences from close relatives for exclusion", args.target_gene)
        excl_recs = get_close_relative_gene_sequences(
            args.organism,
            args.target_gene,
            n=args.max_relative_seqs,
            email=args.email,
            api_key=args.api_key,
        )
    else:
        log.info("Fetching close relatives for exclusion screen")
        excl_recs = get_close_relatives(
            args.organism,
            n=args.max_relative_seqs,
            email=args.email,
            api_key=args.api_key,
        )

    # Auto-cap exclusion sequences for large genomes to keep runtime reasonable
    # If total exclusion sequence size > 10MB, limit to 5 sequences
    total_excl_bp = sum(len(r.seq) for r in excl_recs)
    if total_excl_bp > 10_000_000 and len(excl_recs) > 5:
        log.warning(
            "Large exclusion set detected (%d sequences, %.1fMB total). "
            "Capping to 5 sequences for performance. Use --max-relative-seqs to override.",
            len(excl_recs), total_excl_bp / 1_000_000
        )
        excl_recs = excl_recs[:5]

    # Fetch any extra exclusion accessions explicitly specified — added AFTER cap so they are always included
    if args.extra_exclusion:
        log.info("Fetching %d extra exclusion accession(s): %s", len(args.extra_exclusion), args.extra_exclusion)
        extra_recs = fetch_sequences(args.extra_exclusion, email=args.email, api_key=args.api_key)
        if extra_recs:
            excl_recs.extend(extra_recs)
            log.info("Added %d extra exclusion sequence(s)", len(extra_recs))
        else:
            log.warning("Could not fetch extra exclusion accessions: %s", args.extra_exclusion)

    log.info("Loaded %d exclusion sequence(s)", len(excl_recs))
    return target_recs, excl_recs


def stage_roi(args, target_recs, excl_recs) -> List[ROICandidate]:
    """Stage 2: Find unique ROIs by k-mer subtraction."""
    log = logging.getLogger("stage.roi")

    # Auto-adjust window/step for short sequences (e.g. gene-targeted mode)
    roi_window = args.roi_window
    roi_step   = args.roi_step
    if target_recs:
        median_len = sorted(len(r.seq) for r in target_recs)[len(target_recs) // 2]
        if median_len < roi_window:
            roi_window = max(100, median_len - 20)
            roi_step   = max(10, roi_window // 5)
            log.info("Short sequences detected (median=%dbp) — adjusting window=%d step=%d",
                     median_len, roi_window, roi_step)

    log.info("Scanning for unique ROIs (window=%d step=%d min_uniqueness=%.2f)",
             roi_window, roi_step, args.min_uniqueness)

    rois = find_roi_candidates(
        target_recs,
        excl_recs,
        window_size=roi_window,
        step=roi_step,
        min_uniqueness=args.min_uniqueness,
        top_n=args.top_rois,
    )

    if not rois:
        log.warning("No ROIs found at uniqueness=%.2f. Retrying with 0.60…", args.min_uniqueness)
        rois = find_roi_candidates(
            target_recs,
            excl_recs,
            window_size=roi_window,
            step=roi_step,
            min_uniqueness=0.60,
            top_n=args.top_rois,
        )

    if not rois:
        log.error("Still no ROIs found. Try --min-uniqueness 0.5 or check your sequences.")
        sys.exit(1)

    log.info("Selected %d ROI candidate(s)", len(rois))
    for roi in rois:
        log.info("  %s", roi)

    return rois


def stage_design(args, rois) -> List[PrimerSet]:
    """Stage 3: Design primers/probes with Primer3."""
    log = logging.getLogger("stage.design")

    # Apply CLI overrides to config
    config.PRIMER_MIN_SIZE = args.primer_min_size
    config.PRIMER_OPT_SIZE = args.primer_opt_size
    config.PRIMER_MAX_SIZE = args.primer_max_size
    config.PRIMER_MIN_TM   = args.primer_min_tm
    config.PRIMER_OPT_TM   = args.primer_opt_tm
    config.PRIMER_MAX_TM   = args.primer_max_tm

    all_sets: List[PrimerSet] = []

    for idx, roi in enumerate(rois, 1):
        log.info("Designing primers for ROI %d/%d: %s:%d-%d", idx, len(rois), roi.record_id, roi.start, roi.end)
        sets = design_primers(roi, n_pairs=5)
        log.info("  → %d primer sets designed", len(sets))
        all_sets.extend(sets)

    if not all_sets:
        log.error("Primer3 returned no primer sets for any ROI. "
                  "Try relaxing constraints (--primer-max-size, --primer-min-tm, etc.)")
        sys.exit(1)

    # Sort: passing constraint sets first, then by P3 penalty
    all_sets.sort(key=lambda ps: (not ps.passed_constraints, ps.primer3_penalty))
    log.info("Total primer sets: %d  (%d pass constraints)",
             len(all_sets), sum(1 for ps in all_sets if ps.passed_constraints))

    return all_sets


def stage_blast(args, primer_sets) -> List[PrimerSet]:
    """Stage 4: BLAST validation (local BLAST+)."""
    log = logging.getLogger("stage.blast")

    from pipeline.blast_validator import build_local_db
    db_path = build_local_db(
        organism=args.organism,
        email=args.email,
        api_key=args.api_key,
        db_dir=os.path.join(os.path.expanduser("~"), "blast_db"),
        target_gene=args.target_gene if args.target_gene else None,
    )
    if not db_path:
        log.error("Could not build local BLAST db — skipping BLAST.")
        return primer_sets

    validated = validate_all_sets(
        primer_sets,
        target_organism=args.organism,
        max_sets=args.blast_sets,
        db_path=db_path,
    )
    passed = [ps for ps in validated if ps.blast_pass]
    log.info("BLAST-specific sets: %d / %d", len(passed), args.blast_sets)
    return validated


def stage_report(args, primer_sets):
    """Stage 5: Write output files."""
    log = logging.getLogger("stage.report")
    config.OUTPUT_DIR = args.output_dir

    json_path = write_json(primer_sets, args.organism, args.output_dir)
    tsv_path  = write_tsv(primer_sets, args.organism, args.output_dir)
    txt_path  = write_txt(primer_sets, args.organism, args.output_dir)

    log.info("Output files:")
    log.info("  JSON: %s", json_path)
    log.info("  TSV:  %s", tsv_path)
    log.info("  TXT:  %s", txt_path)

    return txt_path


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = _parse_args()
    _setup_logging(args.verbose)
    log  = logging.getLogger("main")

    log.info("=" * 60)
    log.info("  qPCR Pipeline — target: %s", args.organism)
    log.info("=" * 60)

    # Stage 1 — Fetch sequences
    target_recs, excl_recs = stage_fetch(args)

    # Stage 2 — ROI selection
    rois = stage_roi(args, target_recs, excl_recs)

    # Stage 3 — Primer design
    primer_sets = stage_design(args, rois)

    # Stage 4 — BLAST (optional)
    if not args.no_blast:
        primer_sets = stage_blast(args, primer_sets)
    else:
        log.info("BLAST validation skipped (--no-blast)")

    # Stage 5 — Report
    txt_path = stage_report(args, primer_sets)

    # Print the human-readable report to stdout
    print("\n")
    with open(txt_path) as fh:
        print(fh.read())

    # Exit with non-zero if no sets passed
    fully_passed = [
        ps for ps in primer_sets
        if ps.passed_constraints and (args.no_blast or ps.blast_pass)
    ]
    if not fully_passed:
        log.warning("No primer sets passed all filters. Review output and relax constraints if needed.")
        sys.exit(2)

    log.info("Done. Best set: Pair #%d from ROI %s:%d-%d",
             fully_passed[0].pair_index,
             fully_passed[0].roi.record_id,
             fully_passed[0].roi.start,
             fully_passed[0].roi.end)
    sys.exit(0)


if __name__ == "__main__":
    main()