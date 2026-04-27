# ==============================================================================
# Single-cell & MHC-I Immunopeptidomes Peptide Differential Analysis Pipeline
# Description: Automated Wilcoxon testing for scRNA-seq signature scores and expression
#              Limma-based differential abundance analysis for immunopeptidomics.
# ==============================================================================

rm(list=ls()) # Clean environment

library(tidyverse)
library(Seurat)
library(limma)
library(data.table)
library(writexl)
library(purrr) # For elegant list reductions

# ==============================================================================
# 1. Global Parameters 
# ==============================================================================

# --- 1.1 scRNA-seq Signature Score Settings ---
sc_meta_data      <- scRNA_sigscore_dataset@meta.data   # Seurat object for signature score
score_case        <- "Case"                    # Case
score_ctr         <- "Ctr"                     # Control 
ctypes            <- c("XXX", "XXX")           # Select subsets

# --- 1.2 scRNA-seq FindAllMarkers Settings ---
sc_obj            <- scRNA_dataset             # Seurat object for differentially expression
sc_marker_out     <- "/home/Project/Lung_fibrosis_Immupep/Figure_Data/scRNA-seq_avgdiff"

# --- 1.3 Bulk Immunopeptidomics Settings ---
peptide_dir       <- "/home/Project/Lung_fibrosis_Immupep/Immunopeptidomes/Mmu/raw_data/1.peptide"
limma_out         <- "/home/Project/Lung_fibrosis_Immupep/Figure_Data/Peptide_abundance"

# define group
group_case_name   <- "Case"
group_ctr_name    <- "Ctr"

# Immunopeptidomics data set 
samples_ctr       <- c("S_Ctr1_MHCI.csv", "S_Ctr2_MHCI.csv", "S_Ctr3_MHCI.csv")
samples_case      <- c("S_Case1_MHCI.csv", "S_Case2_MHCI.csv", "S_Case3_MHCI.csv")

# Sample IDs
names_ctr         <- paste0(group_ctr_name, "_", seq_along(samples_ctr))
names_case        <- paste0(group_case_name, "_", seq_along(samples_case))


# ==============================================================================
# 2. scRNA-seq Signature Score Differential Analysis (Wilcoxon)
# ==============================================================================
cat("\nCalculating Wilcoxon statistics for signature scores...\n")

# lapply for all subsets
sigscore_diff_list <- lapply(unique(sc_meta_data$celltype), function(ctype) {
  
  temp <- sc_meta_data %>%
    dplyr::filter(celltype == ctype) %>%
    dplyr::select(all_of(c(score_case, score_ctr))) %>%
    drop_na() # remove NA
  
  # Make sure there are enough cells
  if(nrow(temp) < 3) {
    return(data.frame(celltype = ctype, pvalue = NA, log2_fc = NA))
  }
  
  wilcox_result <- wilcox.test(temp[[score_case]], temp[[score_ctr]], paired = TRUE)
  fold_change   <- mean(temp[[score_case]]) - mean(temp[[score_ctr]])
  
  data.frame(
    celltype = ctype, 
    pvalue   = wilcox_result$p.value, 
    log2_fc  = fold_change
  )
})

Sigscore_diff <- bind_rows(sigscore_diff_list)
print(head(Sigscore_diff))


# ==============================================================================
# 3. scRNA-seq Differential Expression (FindAllMarkers)
# ==============================================================================
cat("\nRunning Seurat FindAllMarkers...\n")

# Create output directory if it does not exist
dir.create(dirname(sc_marker_out), recursive = TRUE, showWarnings = FALSE)

use_TS_markers <- FindAllMarkers(
  sc_obj,
  slot              = "scale.data",
  test.use          = "wilcox",
  min.pct           = 0,
  logfc.threshold   = 0,
  min.cells.feature = 1,
  min.cells.group   = 1,
  random.seed       = 2024
)

saveRDS(use_TS_markers, sc_marker_out)


# ==============================================================================
# 4. Bulk Peptide Abundance Differential Analysis (Limma)
# ==============================================================================
cat("\nRunning Limma differential abundance analysis for peptides...\n")

# --- Intensity object ---
process_peptide_intensity <- function(file_name, sample_id, base_dir) {
  file_path <- file.path(base_dir, file_name)
  
  df <- data.table::fread(file_path) %>%
    dplyr::select(Peptide, Intensity) %>%
    dplyr::distinct(Peptide, .keep_all = TRUE) %>% # unique
    mutate(Intensity = replace_na(Intensity, 0)) %>%
    mutate(Intensity = log2(Intensity + 1))        # Log2 trans
  
  colnames(df) <- c("Peptide", sample_id)
  return(df)
}

# --- processing ---
# purrr::map2 
list_ctr  <- purrr::map2(samples_ctr, names_ctr, ~process_peptide_intensity(.x, .y, peptide_dir))
list_case <- purrr::map2(samples_case, names_case, ~process_peptide_intensity(.x, .y, peptide_dir))

# combine list
all_samples_list <- c(list_ctr, list_case)
intensity_exp_df <- purrr::reduce(all_samples_list, full_join, by = "Peptide") %>%
  replace(is.na(.), 0)

# Convert to expression matrix format
intensity_exp_matrix <- as.matrix(intensity_exp_df[, -1])
rownames(intensity_exp_matrix) <- intensity_exp_df$Peptide

# --- 4.3 Dynamically construct Limma design matrix ---
# Generate factor dynamically based on input sample counts, replacing rep("Ctr", 3)
group_list <- factor(c(
  rep(group_ctr_name, length(samples_ctr)), 
  rep(group_case_name, length(samples_case))
), levels = c(group_ctr_name, group_case_name))
design <- model.matrix(~ group_list) 
rownames(design) <- colnames(intensity_exp_matrix)

# --- 4.4 Differential calculation ---
fit <- lmFit(intensity_exp_matrix, design)
fit <- eBayes(fit, trend = TRUE
result_limma <- topTable(fit, coef = 2, n = Inf) 

# Output
final_result_limma <- result_limma %>% 
  mutate(Peptide = rownames(.)) %>%
  dplyr::select(Peptide, everything())

dir.create(dirname(limma_out), recursive = TRUE, showWarnings = FALSE)
writexl::write_xlsx(final_result_limma, limma_out)

cat("\n==================================================")
cat("\nPipeline Execution Finished Successfully!")
cat("\nLimma results saved to:", limma_out, "\n")
cat("==================================================\n")
