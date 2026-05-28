# qPCR Primer/Probe Design Pipeline

Automated pipeline for designing and validating qPCR primer/probe sets for any target organism, with multiplex compatibility checking.

## Features

- **ROI Selection** — k-mer subtraction against close relatives to find unique genomic regions
- **Primer3 Design** — forward primer, reverse primer, and hydrolysis probe per ROI
- **BLAST Validation** — local BLAST with hit specificity parsing and species-level evaluation
- **Gene-targeted mode** — target a specific gene (e.g. cpn60, rpoB, vly) instead of whole genome
- **Custom FASTA input** — supply your own target sequence directly, bypassing NCBI fetch
- **Extra exclusion accessions** — add specific strains to the exclusion set
- **Relaxable constraints** — Tm, GC, primer size all adjustable via CLI for difficult organisms
- **Multiplex checker** — cross-dimer matrix, Tm spread check, fluorophore channel assignment
- **Three output formats** — JSON (full detail), TSV (spreadsheet), TXT (human-readable)

## Pipeline Stages

```
1. Fetch        — NCBI or local FASTA for target + close relatives
2. ROI Select   — k-mer uniqueness scan (window sliding, deduplication)
3. Primer3      — design FWD/REV/probe for each ROI candidate
4. BLAST        — local blastn, parse top hits, flag non-target matches
5. Report       — JSON + TSV + TXT output
```

## Install

### Quick Setup

```bash
bash setup.sh
```

Handles everything automatically — detects your OS, installs BLAST+, installs Python dependencies, creates required directories.

---

### Manual Setup

#### 1. BLAST+ (system dependency — required)

**Mac:**
```bash
# Option A — conda (recommended)
conda install -c bioconda -c conda-forge blast --override-channels -y

# Option B — Homebrew
brew install blast
```

**Linux / WSL (Windows):**
```bash
sudo apt-get install -y ncbi-blast+
```

Verify:
```bash
blastn -version   # should print: blastn: 2.x.x+
```

#### 2. Python dependencies

```bash
pip install -r requirements.txt
```

#### 3. NCBI credentials

- **Email** — any valid email, required by NCBI for E-utilities
- **API key** — free, optional but recommended (10x rate limit). Get at [ncbi.nlm.nih.gov/account](https://ncbi.nlm.nih.gov/account)

#### 4. Local BLAST database

Auto-built on first run per organism into `~/blast_db/`. Reused on subsequent runs. Delete to force rebuild:
```bash
rm -f ~/blast_db/<organism>_db*
```

---

## Usage

### Standard run (NCBI fetch)
```bash
python run_pipeline.py \
    --organism "Mycolicibacterium fortuitum" \
    --email your@email.com \
    --api-key YOUR_NCBI_API_KEY
```

### Gene-targeted mode
For organisms where whole-genome specificity is insufficient (e.g. Gardnerella).
Fetches only sequences for the specified gene and uses relative gene sequences for exclusion.
```bash
python run_pipeline.py \
    --organism "Gardnerella vaginalis" \
    --target-gene cpn60 \
    --email your@email.com \
    --api-key YOUR_NCBI_API_KEY
```

### Custom FASTA input
Supply your own target sequence directly — useful for specific reference sequences (e.g. vly gene):
```bash
python run_pipeline.py \
    --organism "Gardnerella vaginalis" \
    --target-fasta vly_reference.fasta \
    --email your@email.com \
    --api-key YOUR_NCBI_API_KEY
```

### Extra exclusion accessions
Force specific strains into the exclusion set regardless of size cap — useful when a known cross-reactive species needs to be explicitly excluded:
```bash
python run_pipeline.py \
    --organism "Rickettsia parkeri" \
    --extra-exclusion CP013133.1 CP015012.1 NC_003103.1 CP001612.1 \
    --email your@email.com \
    --api-key YOUR_NCBI_API_KEY
```

### Relaxing constraints for difficult organisms
Some organisms (e.g. Rickettsia) require relaxed Tm or size constraints to get Primer3 to return results:
```bash
python run_pipeline.py \
    --organism "Rickettsia parkeri" \
    --email your@email.com \
    --api-key YOUR_NCBI_API_KEY \
    --primer-min-tm 57 \
    --primer-max-tm 64 \
    --primer-min-size 18 \
    --primer-max-size 28 \
    --top-rois 10
```

### Design-only (skip BLAST)
```bash
python run_pipeline.py \
    --organism "Cyclospora cayetanensis" \
    --email your@email.com \
    --no-blast
```

---

## Multiplex Checker

`multiplex_check.py` takes JSON outputs from individual pipeline runs and checks compatibility for multiplexing multiple targets in a single reaction.

**What it checks:**
- Cross-dimer matrix — all pairwise heterodimer ΔG between every primer and probe across all organisms
- Tm spread — all primer Tms must be within 5°C
- Fluorophore channel assignment

**Separate thresholds for primers vs probes** — probes run at lower concentration so a more lenient threshold (-9.0 kcal/mol) is used for probe-probe pairs vs primers (-6.0 kcal/mol).

**Automatic combination search** — tries all combinations of top-N sets per organism and finds the best combination by fewest cross-dimer failures, then smallest Tm spread.

### Usage
```bash
python multiplex_check.py \
    --panel \
        "FAM:output/mycolicibacterium_fortuitum_results.json" \
        "SUN:output/mycobacteroides_chelonae_results.json" \
        "CY5:output/cyclospora_cayetanensis_results.json" \
        "TXRED:output/homo_sapiens_results.json" \
    --output output/multiplex_report.txt \
    --top-n 3
```

### Multiplex parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--top-n` | 3 | Top N sets per organism to try |
| `--dimer-threshold` | -6.0 | Primer heterodimer ΔG threshold (kcal/mol) |
| `--probe-dimer-threshold` | -9.0 | Probe heterodimer ΔG threshold (kcal/mol) |
| `--tm-spread` | 5.0 | Max allowed Tm spread (°C) |

---

## Key Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--primer-min-size` | 19 | Min primer length |
| `--primer-max-size` | 26 | Max primer length |
| `--primer-min-tm` | 59.0 | Min primer Tm (°C) |
| `--primer-max-tm` | 62.0 | Max primer Tm (°C) |
| `--roi-window` | 500 | ROI window size (bp) |
| `--min-uniqueness` | 0.80 | Min k-mer uniqueness (0–1) |
| `--blast-sets` | 3 | Number of primer sets to BLAST validate |
| `--top-rois` | 5 | Number of top ROIs to design primers for |
| `--max-relative-seqs` | 10 | Max close relative sequences for exclusion |
| `--extra-exclusion` | — | Extra NCBI accessions to add to exclusion set |
| `--target-gene` | — | Target a specific gene instead of whole genome |
| `--target-fasta` | — | Use a local FASTA as target (skips NCBI fetch) |
| `--no-blast` | — | Skip BLAST validation (design-only) |

## Constraints (from config.py)

- Primer length: 19–26 nt (adjustable via CLI)
- Primer Tm: 59–62°C (adjustable via CLI)
- Probe length: 22–29 nt
- Probe Tm: 5–10°C above mean primer Tm
- Amplicon: 70–200 bp
- Self-dimer ΔG: ≥ -6.0 kcal/mol

## Output Files (in `./output/`)

- `<organism>_results.json` — full pipeline output, machine-readable
- `<organism>_results.tsv` — one row per primer set, spreadsheet-ready
- `<organism>_report.txt` — human-readable summary with BLAST hit tables
- `multiplex_report.txt` — multiplex compatibility report (from multiplex_check.py)

## Validated Organisms

| Organism | Mode | Result |
|----------|------|--------|
| Mycolicibacterium fortuitum | Whole genome | ✅ Specific |
| Mycobacteroides chelonae | Whole genome | ✅ Specific |
| Cyclospora cayetanensis | Whole genome | ✅ Specific |
| Gardnerella vaginalis | vly gene (--target-fasta) | ✅ Specific (79% sensitivity — vly-carrying strains only) |
| Rickettsia parkeri | Whole genome + extra-exclusion | ✅ Specific (insect endosymbiont hits only, 4-5 mismatches, no clinical relevance) |
| Rickettsia rickettsii | Whole genome | ❌ Unique regions are repetitive elements — unprimable. A1G_04230 is 99% identical to R. lanei |

## Exit Codes

- `0` — at least one primer set passed all constraints + BLAST
- `2` — pipeline ran but no sets fully passed (relax constraints)
- `1` — fatal error (no sequences found, etc.)
