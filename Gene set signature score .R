# ==============================================================================
# Signature Scoring Pipeline
# Description: Calculates signature scores for defined gene sets across single-cell
#              populations and exports the results as a clean CSV file.
# ==============================================================================

rm(list=ls()) # Clean environment

library(tidyverse)
library(Seurat)
library(pagoda2)

# --- 1. Define Global Parameters ---

# Input & Output Paths
scRNA_rds_path <- "/home/Project/Lung_fibrosis_Immupep/scRNA-seq/scRNA_data.rds"
output_dir <- "/home/Project/Lung_fibrosis_Immupep/Figure_Data/Signature_score"
output_filename <- "G_Signature_Score.pdf"


# Define Gene Sets (Dynamically expandable)
# Ensure PF_unique_protein and Con_unique_protein are loaded before running
gene_sets <- list(
  "GeneSet1" = GeneSet1$Symbol,
  "GeneSet2" = GeneSet2$Symbol
)

# Select Cell Types (Optional: Set to NULL if you want to output all cells)
ctypes <- c("XXX", "XXX")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- 2. Load and Prepare Data ---
cat("[1/4] Loading scRNA-seq object...\n")
sc_obj <- readRDS(scRNA_rds_path)

# Extract count matrix and metadata
countMatrix <- t(sc_obj@assays$RNA@counts)
meta <- sc_obj@meta.data
meta$Barcode <- rownames(meta)

# --- 3. Calculate Signature Scores ---
cat("[2/4] Calculating signature scores...\n")

# Iterate over the gene sets and calculate scores
score_list <- lapply(names(gene_sets), function(gs_name) {
  scores <- score.cells.puram(countMatrix, gene_sets[[gs_name]])
  data.frame(
    Barcode = names(scores),
    Score = scores,
    Signature = gs_name
  )
})

# Combine and convert to Wide-format (Rows: Cells, Columns: Signature Scores)
cat("[3/4] Formatting data...\n")
score_df_long <- bind_rows(score_list)
score_df_wide <- score_df_long %>% 
  pivot_wider(names_from = Signature, values_from = Score)

# Merge with celltype metadata
final_results <- meta %>%
  dplyr::select(Barcode, celltype) %>%
  inner_join(score_df_wide, by = "Barcode")

# Filter for target cell types
if (!is.null(target_celltypes)) {
  final_results <- final_results %>% filter(celltype %in% target_celltypes)
}

# --- 4. Export Results ---
cat("[4/4] Exporting results to CSV...\n")
output_path <- file.path(output_dir, output_csv_name)

write.csv(final_results, file = output_path, row.names = FALSE, quote = FALSE)

cat("==================================================\n")
cat("Done! Results successfully saved to:\n")
cat(output_path, "\n")
cat("==================================================\n")
