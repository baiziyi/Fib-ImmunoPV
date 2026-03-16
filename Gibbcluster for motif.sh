#!/bin/bash
# ==============================================================================
# GibbsCluster Automated Pipeline
# Description: Batch unsupervised clustering of immunopeptidomics data.
#              Optimized to prevent data loss of 8-mer peptides.
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -e

# --- 1. Define Global Parameters ---
INPUT_DIR="/home/Project/Lung_fibrosis_Immupep/Figure_Data/filter_Gibbcluster"   # Directory containing input peptide lists (.txt)
OUTPUT_DIR="/home/Project/Lung_fibrosis_Immupep/Figure_Data/Gibbs"        # Main directory for all GibbsCluster outputs
GIBBS_BIN="/home/resource/software/gibbscluster-2.0/gibbscluster"

# Core GibbsCluster Settings
MOTIF_LENGTH=9         # (-l) Base motif length
CLUSTERS="1-5"         # (-g) Number of clusters to test
MAX_INSERTION=2        # (-I) Max insertion length (rescues 10-11mers)
MAX_DELETION=1         # (-D) Max deletion length (rescues 8-mers!)
TRASH_CLUSTER=1        # (-t) Enable trash cluster to remove outliers (highly recommended)

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=================================================="
echo "Starting GibbsCluster Pipeline..."
echo "Input directory:  ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Motif Length:     ${MOTIF_LENGTH} (with Insertions: ${MAX_INSERTION}, Deletions: ${MAX_DELETION})"
echo "=================================================="

# --- 2. Iterate through all peptide files ---
for PEP_FILE in "${INPUT_DIR}"/*.txt; do
    
    # Safely check if files exist
    if [[ ! -f "${PEP_FILE}" ]]; then
        echo "[WARNING] No .txt peptide files found in ${INPUT_DIR}. Exiting."
        break
    fi

    # Extract base sample name safely
    FILENAME=$(basename "${PEP_FILE}")
    SAMPLE="${FILENAME%.txt}"
    
    echo -e "\n>>> Clustering Sample: ${SAMPLE}"

    # Create sample-specific output directory to avoid clutter
    # GibbsCluster generates many files, so isolating them per sample is crucial
    SAMPLE_OUT_DIR="${OUTPUT_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_OUT_DIR}"
    
    # Switch to the output directory before running (GibbsCluster tends to write locally)
    pushd "${SAMPLE_OUT_DIR}" > /dev/null

    # --- Step 3: Execute GibbsCluster ---
    # We use absolute path for the input file since we changed directories
    ABS_PEP_FILE=$(readlink -f "../../${PEP_FILE}")

    # Note: Added -l, -g, -I, -D, -t flags to optimize for variable-length MHC-I peptides
    ${GIBBS_BIN} \
        -f "${ABS_PEP_FILE}" \
        -P "${SAMPLE}" \
        -l ${MOTIF_LENGTH} \
        -g ${CLUSTERS} \
        -I ${MAX_INSERTION} \
        -D ${MAX_DELETION} \
        -t ${TRASH_CLUSTER} \
        > "${SAMPLE}.stdout.log"

    popd > /dev/null

    echo "    Successfully completed: ${SAMPLE}"
    echo "    Results saved to: ${SAMPLE_OUT_DIR}"
done

echo -e "\n=================================================="
echo "All clustering jobs finished successfully!"
echo "Please check the KLD scores in the output logs to select the optimal cluster number."
echo "=================================================="