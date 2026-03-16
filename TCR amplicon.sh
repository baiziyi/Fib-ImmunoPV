#!/bin/bash
# ==============================================================================
# MiXCR TCR Amplicon Processing Pipeline
# Description: Batch processing of paired-end TCR-seq FASTQ files for mouse (mmu).
# Usage: bash run_mixcr_pipeline.sh
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -e

# --- 1. Define Global Parameters  ---
INPUT_DIR="."                     # Directory containing raw fastq.gz files
OUTPUT_DIR="../TRB_mixcr4"        # Main output directory
SPECIES="mmu"                     # Species (mmu for mouse, hsa for human)
MEMORY="40g"                      # Java heap memory limit
THREADS=8                         # Number of threads (optional, if you want to limit CPU usage)

# Create main output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "=================================================="
echo "Starting MiXCR Pipeline..."
echo "Input directory:  ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Species:          ${SPECIES}"
echo "=================================================="

# --- 2. Iterate through all R1 files ---
# This avoids the fragile `ls | awk` method
for R1_FILE in "${INPUT_DIR}"/*.dedup.R1.fastq.gz; do
    
    # Extract the base sample name using bash string manipulation
    FILENAME=$(basename "${R1_FILE}")
    SAMPLE="${FILENAME%%.dedup.R1.fastq.gz}"
    R2_FILE="${INPUT_DIR}/${SAMPLE}.dedup.R2.fastq.gz"

    # Check if the corresponding R2 file exists
    if [[ ! -f "${R2_FILE}" ]]; then
        echo "[WARNING] Matching R2 file not found for ${SAMPLE}. Skipping..."
        continue
    fi

    echo -e "\n>>> Processing Sample: ${SAMPLE}"

    # Create sample-specific output directory
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_OUT}"
    
    # Define common prefix for all output files of this sample
    PREFIX="${SAMPLE_OUT}/${SAMPLE}"

    # --- Step 1: Alignment ---
    # Matches the exact parameters from your report log
    echo "    [1/3] Aligning reads..."
    mixcr align -f \
        --preset generic-tcr-amplicon \
        --species "${SPECIES}" \
        --rna \
        --rigid-left-alignment-boundary \
        --floating-right-alignment-boundary C \
        --report "${PREFIX}.align.report.txt" \
        --json-report "${PREFIX}.align.report.json" \
        -Xmx"${MEMORY}" \
        "${R1_FILE}" "${R2_FILE}" "${PREFIX}.vdjca"

    # --- Step 2: Assembly ---
    # Incorporates your specific custom flags (-O)
    echo "    [2/3] Assembling clonotypes (CDR3)..."
    mixcr assemble -f \
        --report "${PREFIX}.assemble.report.txt" \
        --json-report "${PREFIX}.assemble.report.json" \
        -OassemblingFeatures="CDR3" \
        -OseparateByJ=true \
        -OseparateByV=true \
        -Xmx"${MEMORY}" \
        "${PREFIX}.vdjca" "${PREFIX}.clns"

    # --- Step 3: Export Clones ---
    # Essential for downstream analysis in R/Python (converts .clns to .txt)
    echo "    [3/3] Exporting clonotype table..."
    mixcr exportClones -f \
        -o -t \
        "${PREFIX}.clns" "${PREFIX}.clones.txt"

    echo "    Successfully completed: ${SAMPLE}"
done

echo -e "\n=================================================="
echo "All samples processed successfully! Data saved to ${OUTPUT_DIR}"
echo "=================================================="