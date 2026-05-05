# qPCR Primer/Probe Design Pipeline

Automated pipeline for designing and validating qPCR primer/probe sets for any target organism.

## Features

- **ROI Selection** — k-mer subtraction against close relatives to find unique genomic regions
- **Primer3 Design** — forward primer, reverse primer, and hydrolysis probe per ROI
- **BLAST Validation** — parallel NCBI BLAST (5 workers) with hit specificity parsing
- **Flexible constraints** — all Tm, GC, length, and amplicon parameters configurable via CLI
- **Three output formats** — JSON (full detail), TSV (spreadsheet), TXT (human-readable)

## Pipeline Stages

```
1. Fetch        — NCBI or local FASTA for target + close relatives
2. ROI Select   — k-mer uniqueness scan (window sliding, deduplication)
3. Primer3      — design FWD/REV/probe for each ROI candidate
4. BLAST        — parallel blastn, parse top 15 hits, flag non-target matches
5. Report       — JSON + TSV + TXT output
```

## Install

### 1. BLAST+ (system dependency — required)

BLAST+ is a compiled binary, not a Python package, so it must be installed separately.

**Mac:**
```bash
# Option A — conda (recommended)
conda install -c bioconda -c conda-forge blast --override-channels -y

# Option B — Homebrew
brew install blast
```

**Linux:**
```bash
conda install -c bioconda blast -y
# or
sudo apt-get install ncbi-blast+
```

Verify:
```bash
blastn -version   # should print: blastn: 2.x.x+
```

### 2. Python dependencies

```bash
pip install -r requirements.txt
```

### 3. NCBI credentials

- **Email** — any valid email, required by NCBI for E-utilities
- **API key** — free, optional but recommended (10x rate limit). Get at [ncbi.nlm.nih.gov/account](https://ncbi.nlm.nih.gov/account)

### 4. Local BLAST database

Auto-built on first run per organism into `~/blast_db/`. Reused on subsequent runs.

## Usage

### Minimal (NCBI fetch)
```bash
python run_pipeline.py \
    --organism "Cyclospora cayetanensis" \
    --email your@email.com
```

### With NCBI API key (10x rate limit)
```bash
python run_pipeline.py \
    --organism "Cyclospora cayetanensis" \
    --email your@email.com \
    --api-key YOUR_NCBI_API_KEY
```

### Relax primer length for difficult organisms
```bash
python run_pipeline.py \
    --organism "Cyclospora cayetanensis" \
    --email your@email.com \
    --primer-max-size 26
```

### Design-only (skip BLAST)
```bash
python run_pipeline.py \
    --organism "Cyclospora cayetanensis" \
    --email your@email.com \
    --no-blast
```

### Pre-downloaded FASTAs (avoids NCBI fetch)
```bash
python run_pipeline.py \
    --organism "Cyclospora cayetanensis" \
    --email your@email.com \
    --target-fasta target.fasta \
    --exclusion-fasta relatives.fasta
```

## Key Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--primer-min-size` | 19 | Min primer length |
| `--primer-max-size` | 26 | Max primer length |
| `--primer-min-tm` | 59.0 | Min primer Tm (°C) |
| `--primer-max-tm` | 62.0 | Max primer Tm (°C) |
| `--roi-window` | 500 | ROI window size (bp) |
| `--min-uniqueness` | 0.80 | Min k-mer uniqueness (0–1) |
| `--blast-workers` | 5 | Parallel BLAST threads |
| `--blast-sets` | 3 | # of primer sets to BLAST |
| `--no-blast` | — | Skip BLAST (design-only) |

## Constraints (from config.py)

- Primer length: 19–26 nt (avoid >22 if possible)
- Primer Tm: 59–62°C
- Probe length: 22–29 nt  
- Probe Tm: 5–10°C above mean primer Tm
- Amplicon: 70–200 bp

## Output Files (in `./output/`)

- `<organism>_results.json` — full pipeline output, machine-readable
- `<organism>_results.tsv` — one row per primer set, spreadsheet-ready
- `<organism>_report.txt` — human-readable summary with BLAST hit tables

## Exit Codes

- `0` — at least one primer set passed all constraints + BLAST
- `2` — pipeline ran but no sets fully passed (relax constraints)
- `1` — fatal error (no sequences found, etc.)
