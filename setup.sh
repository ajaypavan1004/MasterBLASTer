#!/bin/bash
# qPCR Pipeline Setup Script
set -e

echo "=== qPCR Primer/Probe Pipeline Setup ==="

# 1. Check for BLAST+
if command -v blastn &> /dev/null; then
    echo "✓ BLAST+ already installed: $(blastn -version 2>&1 | head -1)"
else
    echo "Installing BLAST+..."
    if command -v conda &> /dev/null; then
        conda install -c bioconda -c conda-forge blast --override-channels -y
    elif command -v brew &> /dev/null; then
        brew install blast
    elif command -v apt-get &> /dev/null; then
        sudo apt-get install -y ncbi-blast+
    else
        echo "ERROR: Could not find conda, brew, or apt-get. Please install BLAST+ manually."
        echo "  Mac:   brew install blast"
        echo "  Linux: sudo apt-get install ncbi-blast+"
        exit 1
    fi
    echo "✓ BLAST+ installed"
fi

# 2. Install Python dependencies
echo "Installing Python dependencies..."
pip install -r requirements.txt -q
echo "✓ Python dependencies installed"

# 3. Create output and db directories
mkdir -p output ~/blast_db
echo "✓ Directories created"

echo ""
echo "=== Setup complete ==="
echo ""
echo "Run the pipeline with:"
echo "  python run_pipeline.py --organism \"Cyclospora cayetanensis\" --email your@email.com"
