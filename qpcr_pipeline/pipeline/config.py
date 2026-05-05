"""
Pipeline configuration — all tunable constraints live here.
"""

# ── Primer constraints ────────────────────────────────────────────────────────
PRIMER_MIN_SIZE     = 19
PRIMER_OPT_SIZE     = 20
PRIMER_MAX_SIZE     = 26          # relaxed to 26 for difficult organisms

PRIMER_MIN_TM       = 59.0
PRIMER_OPT_TM       = 60.0
PRIMER_MAX_TM       = 62.0

PRIMER_MIN_GC       = 40.0
PRIMER_OPT_GC       = 50.0
PRIMER_MAX_GC       = 65.0

# ── Probe constraints ─────────────────────────────────────────────────────────
PROBE_MIN_SIZE      = 22
PROBE_OPT_SIZE      = 25
PROBE_MAX_SIZE      = 29

PROBE_TM_DELTA_MIN  = 5.0        # probe Tm must be this much ABOVE primer Tm
PROBE_TM_DELTA_MAX  = 10.0

PROBE_MIN_GC        = 40.0
PROBE_MAX_GC        = 70.0

# ── Amplicon size ─────────────────────────────────────────────────────────────
AMPLICON_MIN        = 70
AMPLICON_OPT        = 120
AMPLICON_MAX        = 200

# ── BLAST settings ────────────────────────────────────────────────────────────
BLAST_DB            = "refseq_select"
BLAST_PROGRAM       = "blastn"
BLAST_HITLIST       = 15          # top N hits to evaluate
BLAST_WORKERS       = 5           # parallel BLAST threads
BLAST_EVALUE        = 10          # liberal — specificity determined by hit names
BLAST_WORD_SIZE     = 7           # short word size for primers

# ── ROI selection ─────────────────────────────────────────────────────────────
ROI_MIN_LENGTH      = 300         # minimum unique region to feed Primer3
ROI_MAX_WINDOWS     = 20          # max k-mer windows to probe for uniqueness
KMER_K              = 20          # k-mer size for uniqueness screen

# ── NCBI E-utilities ─────────────────────────────────────────────────────────
ENTREZ_EMAIL        = "user@example.com"   # set via CLI --email
ENTREZ_API_KEY      = ""                   # optional; set via CLI --api-key
ENTREZ_TOOL         = "qpcr_pipeline"

# ── Output ────────────────────────────────────────────────────────────────────
OUTPUT_DIR          = "output"
