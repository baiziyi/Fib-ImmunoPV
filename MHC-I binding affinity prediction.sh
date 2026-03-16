#!/bin/bash
# ==============================================================================
# NetMHCpan & MhcVizPipe Automated Pipeline
# Description: Batch processing for MHC-I binding affinity predictions and 
#              automated generation of MhcVizPipe configuration.
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -e

# --- 1. Define Global Parameters  ---
INPUT_DIR="/home/Project/Lung_fibrosis_Immupep/Figure_Data/Immunopeptidomes/Hsa/raw_data/1.peptide"      # Directory containing peptide list files (.pep or .txt)
OUTPUT_DIR="/home/Project/Lung_fibrosis_Immupep/Figure_Data/NetMHCpan"   # Main directory for prediction results
MVP_CONFIG_FILE="./mvp_workspace/config.ini" # Path to generate MhcVizPipe config

# Core NetMHCpan settings
# Note: Defaulting to C57BL/6 murine alleles. Change to 'HLA-A02:01...' for human data.
ALLELES="H-2-Kb,H-2-Db"           
PEPTIDE_LENGTHS="8,9,10,11,12,13,14,15,16" # Matches your "class I max length = 16"
THRESHOLD_BINDER=2.0              # %Rank threshold for weak binders
THRESHOLD_STRONG=0.5              # %Rank threshold for strong binders

# Create necessary directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname "${MVP_CONFIG_FILE}")"

echo "=================================================="
echo "Starting NetMHCpan Prediction Pipeline..."
echo "Input directory:  ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Target Alleles:   ${ALLELES}"
echo "Max Length:       16 aa"
echo "=================================================="

# --- 2. Iterate through all peptide files ---
for PEP_FILE in "${INPUT_DIR}"/*.txt; do
    
    # Check if files exist to prevent globbing errors
    if [[ ! -f "${PEP_FILE}" ]]; then
        echo "[WARNING] No .txt peptide files found in ${INPUT_DIR}."
        break
    fi

    # Extract base sample name
    FILENAME=$(basename "${PEP_FILE}")
    SAMPLE="${FILENAME%.txt}"
    
    echo -e "\n>>> Predicting affinities for Sample: ${SAMPLE}"
    
    # Define outputs
    SAMPLE_OUT_TXT="${OUTPUT_DIR}/${SAMPLE}_netMHCpan.txt"
    SAMPLE_OUT_XLS="${OUTPUT_DIR}/${SAMPLE}_netMHCpan.xls"

    # Run NetMHCpan command
    # -p indicates input is a peptide list (not FASTA)
    # -BA includes Binding Affinity predictions
    # -xls saves a highly readable table for downstream R/Python parsing
    netMHCpan -p "${PEP_FILE}" \
        -a "${ALLELES}" \
        -l "${PEPTIDE_LENGTHS}" \
        -th "${THRESHOLD_BINDER}" \
        -sth "${THRESHOLD_STRONG}" \
        -BA \
        -xls \
        -xlsfile "${SAMPLE_OUT_XLS}" > "${SAMPLE_OUT_TXT}"

    echo "    Successfully completed: ${SAMPLE}"
done

echo -e "\n=================================================="
echo "Generating MhcVizPipe Configuration..."
echo "=================================================="

# --- 3. Auto-generate MhcVizPipe config.ini ---
# This matches your exact requested parameters and writes them to the config file
cat << EOF > "${MVP_CONFIG_FILE}"
[DIRECTORIES]
NetMHCpan path = AUTO
GibbsCluster path = AUTO
temp directory = tmpmhcvizpipe

[ANALYSIS]
motifs = kullback-leibler
hobohm clustering = yes
clustering threshold = 0.63
weight on prior = 200
max threads = -1
class I max length = 16

[SERVER]
HOSTNAME = 0.0.0.0
PORT = 8080
EOF

echo "MhcVizPipe config saved to: ${MVP_CONFIG_FILE}"
echo "MHC-I binding affinity predictions pipeline execution finished gracefully!"
