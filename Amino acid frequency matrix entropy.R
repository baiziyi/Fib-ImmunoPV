# ==============================================================================
# Immunopeptidome Motif Distance and Entropy Pipeline
# ==============================================================================

rm(list=ls())

library(tidyverse)
library(data.table)
library(entropy)

# ==============================================================================
# 1. Global Parameters
# ==============================================================================

group_ctrl <- "Ctr"
group_case <- "Case"

# Directories
dir_excel   <- "/home/Project/Lung_fibrosis_Immupep/Immunopeptidomes/Hsa/raw_data/1.peptide"
dir_pep_out <- "/home/Project/Lung_fibrosis_Immupep/Figure_Data/scRNA-seq_avgdiff"
dir_gibbs   <- "/home/Project/Lung_fibrosis_Immupep/Figure_Data/Gibbs"
dir_csv_out <- "/home/Project/Lung_fibrosis_Immupep/Figure_Data/aaFreq_matrix"

dir.create(dir_pep_out, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_csv_out, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Dynamic Peptide Extraction (Excel to TXT)
# ==============================================================================

regex_pattern <- paste0("^(", group_ctrl, "|", group_case, ").*\\.xlsx$")
excel_files   <- list.files(dir_excel, pattern = regex_pattern, full.names = TRUE)

for (file in excel_files) {
  sample_name <- tools::file_path_sans_ext(basename(file))
  df <- readxl::read_xlsx(file)
  
  out_path <- file.path(dir_pep_out, paste0(sample_name, "_pep.txt"))
  write.table(unique(df[, "Peptide"]), out_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# ==============================================================================
# 3. Motif Matrix Loading & Normalization
# ==============================================================================

sample_dirs <- list.dirs(dir_gibbs, recursive = FALSE, full.names = FALSE)
target_samples <- sample_dirs[grepl(paste0("^(", group_ctrl, "|", group_case, ")"), sample_dirs)]

motif_list <- list()

for (sample in target_samples) {
  file_path <- file.path(dir_gibbs, sample, "logos/gibbs_logos_1of1_freq.mat")
  if (file.exists(file_path)) {
    mat <- fread(file_path)
    motif_list[[sample]] <- as.matrix(mat[, 3:22]) 
  }
}

motif_norm <- lapply(motif_list, function(mat) {
  t(apply(mat, 1, function(row) row / sum(row)))
})

# ==============================================================================
# 4. Pairwise Distance Calculation & Export
# ==============================================================================

calc_distance <- function(mat1, mat2) {
  dist_matrix <- mapply(function(r1, r2) sqrt(sum((r1 - r2)^2)), 
                        split(mat1, row(mat1)), split(mat2, row(mat2)))
  return(mean(dist_matrix))
}

n_samples <- length(motif_norm)
sample_names <- names(motif_norm)
dist_mat <- matrix(0, n_samples, n_samples, dimnames = list(sample_names, sample_names))

for (i in 1:n_samples) {
  for (j in i:n_samples) {
    d <- calc_distance(motif_norm[[i]], motif_norm[[j]])
    dist_mat[i, j] <- d
    dist_mat[j, i] <- d
  }
}

# Export distance matrix
dist_out_path <- file.path(dir_csv_out, "Motif_Distance_Matrix.csv")
write.csv(dist_mat, dist_out_path, quote = FALSE)

# ==============================================================================
# 5. Positional Entropy Calculation & Export
# ==============================================================================

entropy_res <- lapply(motif_norm, function(mat) apply(mat, 1, entropy.empirical))
entropy_df <- as.data.frame(do.call(cbind, entropy_res))
entropy_df$Position <- paste0("P", 1:nrow(entropy_df))

# Relocate 'Position' to be the first column
entropy_df <- entropy_df %>% relocate(Position)

# Export wide format
entropy_wide_path <- file.path(dir_csv_out, "Positional_Entropy_Wide.csv")
write.csv(entropy_df, entropy_wide_path, row.names = FALSE, quote = FALSE)

# Reshape to long format and export
long_df <- pivot_longer(entropy_df, cols = -Position, names_to = "Sample", values_to = "Entropy") %>%
  mutate(Group = ifelse(grepl(paste0("^", group_ctrl), Sample), group_ctrl, group_case))

entropy_long_path <- file.path(dir_csv_out, "Positional_Entropy_Long.csv")
write.csv(long_df, entropy_long_path, row.names = FALSE, quote = FALSE)
