# ==============================================================================
# Fib-SCORE Pipeline: Multi-Omics Immunopeptidome Prioritization
# Description: Integrates stringent quality control, presentation stability, 
#              allelic binding potential, and targeting specificity to 
#              prioritize actionable therapeutic peptides in fibrosis.
# ==============================================================================

rm(list=ls())

library(tidyverse)
library(data.table)
library(writexl)

# ==============================================================================
# 1. Global Parameters
# ==============================================================================

group_ctr  <- "Ctr"
group_case <- "Case"

# Criterion (i): Peptide identification confidence & (ii) Predicted binding affinity
ms2_cor_thresh   <- 0.4
caa_thresh       <- 60
delta_rt_thresh  <- 30
aff_nm_thresh    <- 1000

# Criterion (iiii): Sample presentation stability
min_case_count   <- 2    # Detected in more than half of samples

# Target Cell Types for Target Score (Pathological cell populations)
target_celltypes <- c("IM", "SAM_int", "SAM", 
                      "Mesothelial", "FB_Matrix", "MyoFB_Alv", "Endo")

# Output path
out_score_file   <- "/home/Project/Lung_fibrosis_Immupep/Figure_Data/Fib_SCORE/Peptide_score_Integrated.xlsx"
dir.create(dirname(out_score_file), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Input Data Placeholders (Replace with your actual loaded objects)
# ==============================================================================
# Input 1: Mmu_abundance_df (Combined dataframe of all samples)
# Must contain columns: Sample, Peptide, Accession, Mmu_Symbol, Delta RT, 
# MS2 Correlation, CAA (%), Aff_nM, Aff_Rank

# Input 2: Hsa_peptide_geneset (Human encoded the source proteins for the peptides list)
# e.g., Hsa_peptide_geneset <- data.frame(Symbol = c("..."))

# Input 3: Human_mouse_gene_maping (Ortholog mapping)
# e.g., Human_mouse_gene_maping <- data.frame(Mmu_Symbol = c("..."), Symbol = c("..."))

# Input 4: Mmu_result_limma (Limma differential abundance results)
# Must contain: Peptide, logFC, adj.P.Val

# Input 5: sc_markers (scRNA_data.markers, FindAllMarkers with scale.data)
# Must contain: Mmu_Symbol, cluster, avg_diff

# ==============================================================================
# 3. Filtering Step 1: Identification Confidence & Presentation Stability
# ==============================================================================
cat("[1/5] Applying Peptide Identification Confidence & Stability filters...\n")

# Assign group based on sample name
Mmu_abundance_df <- Mmu_abundance_df %>%
  mutate(Group = ifelse(grepl(paste0("^", group_ctr), Sample), group_ctr, group_case))

# Apply mass spec (identification confidence) & affinity filters
filtered_abundance <- Mmu_abundance_df %>%
  filter(abs(`Delta RT`) < delta_rt_thresh & 
           `MS2 Correlation` > ms2_cor_thresh & 
           `CAA (%)` > caa_thresh & 
           Aff_nM < aff_nm_thresh)

# Calculate occurrence count in Case vs Ctr (Presentation stability)
pep_counts <- filtered_abundance %>%
  group_by(Peptide, Group) %>%
  summarise(n_samples = n_distinct(Sample), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = n_samples, values_fill = 0)

# Ensure columns exist even if one group is completely missing
if(!group_ctr %in% names(pep_counts)) pep_counts[[group_ctr]] <- 0
if(!group_case %in% names(pep_counts)) pep_counts[[group_case]] <- 0

# Base candidate filter: Must be present in >= half of case samples
candidates_base <- pep_counts %>%
  filter(!!sym(group_case) >= min_case_count)

# Extract unique mapping (Peptide to Mmu_Symbol)
pep_to_gene <- filtered_abundance %>%
  dplyr::select(Peptide, Mmu_Symbol) %>%
  distinct() %>%
  drop_na(Mmu_Symbol)

# ==============================================================================
# 4. Filtering Step 2: Source Protein Conservation
# ==============================================================================
cat("[2/5] Mapping orthologs to ensure source protein conservation...\n")

# Map human target genes to mouse symbols (mapped homologs in humans)
target_mmu_genes <- Hsa_peptide_geneset %>%
  inner_join(Human_mouse_gene_maping, by = "Symbol") %>%
  dplyr::select(Mmu_Symbol) %>%
  distinct()

# Filter candidates by human ortholog overlap
candidates_overlap <- candidates_base %>%
  inner_join(pep_to_gene, by = "Peptide") %>%
  inner_join(target_mmu_genes, by = "Mmu_Symbol")

# ==============================================================================
# 5. Fib-SCORE Calculation 
# ==============================================================================

# --- 5.1 Specificity Score (Peptide expression abundance differential) ---
cat("[3/5] Calculating Specificity (Abundance) Score...\n")
# Prioritizes peptides that are both disease-specific and stably presented

abundance_score_df <- candidates_overlap %>%
  left_join(Mmu_result_limma, by = "Peptide") %>%
  mutate(
    Abundance_score = case_when(
      !!sym(group_ctr) == 0 ~ 1,                         # Unique to Case
      logFC > 0 ~ 1 - adj.P.Val,                         # Upregulated
      logFC < 0 ~ -1 * (1 - adj.P.Val),                  # Downregulated
      TRUE ~ 0
    )
  ) %>%
  dplyr::select(Peptide, Mmu_Symbol, Abundance_score) %>%
  replace_na(list(Abundance_score = 0))

# --- 5.2 Target Score (Potential targeting of pathological cell populations) ---
cat("[4/5] Calculating Target Score (scRNA-seq Fisher's Test)...\n")
# Prioritizes peptides targeting homologous pathological cell types

unique_candidate_genes <- unique(abundance_score_df$Mmu_Symbol)

# Define Target vs Other clusters
all_clusters   <- unique(sc_markers$cluster)
other_clusters <- setdiff(all_clusters, target_celltypes)

# Count enriched (avg_diff > 0) vs unenriched (<= 0) across Target vs Other
fisher_df <- sc_markers %>%
  filter(Mmu_Symbol %in% unique_candidate_genes) %>%
  mutate(
    is_target_cluster = cluster %in% target_celltypes,
    is_enriched       = avg_diff > 0
  ) %>%
  group_by(Mmu_Symbol) %>%
  summarise(
    en_target   = sum(is_target_cluster & is_enriched, na.rm = TRUE),
    en_other    = sum(!is_target_cluster & is_enriched, na.rm = TRUE),
    unen_target = sum(is_target_cluster & !is_enriched, na.rm = TRUE),
    unen_other  = sum(!is_target_cluster & !is_enriched, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate Target score via rowwise Fisher test
fisher_df <- fisher_df %>%
  rowwise() %>%
  mutate(
    Fisher_pvalue = fisher.test(matrix(c(en_target, unen_target, 
                                         en_other, unen_other), nrow = 2))$p.value,
    Target_score  = 1 - Fisher_pvalue
  ) %>%
  ungroup() %>%
  dplyr::select(Mmu_Symbol, Target_score)

# --- 5.3 Affinity Score (Predicted peptide-MHC-I binding affinity) ---
cat("[5/5] Calculating Affinity Score...\n")
# Prioritizes peptides with high allelic binding potential

affinity_score_df <- filtered_abundance %>%
  filter(Peptide %in% abundance_score_df$Peptide) %>%
  mutate(Affinity_score = (100 - Aff_Rank) / 100) %>%
  group_by(Peptide) %>%
  slice_max(Affinity_score, n = 1, with_ties = FALSE) %>% # Get best score per peptide
  ungroup() %>%
  dplyr::select(Peptide, Affinity_score)

# ==============================================================================
# 6. Final Fib-SCORE Integration and Export
# ==============================================================================

final_scores <- abundance_score_df %>%
  left_join(fisher_df, by = "Mmu_Symbol") %>%
  left_join(affinity_score_df, by = "Peptide") %>%
  replace_na(list(Target_score = 0, Affinity_score = 0)) %>%
  mutate(Final_Score = Abundance_score + Target_score + Affinity_score) %>%
  arrange(desc(Final_Score)) %>%
  group_by(Mmu_Symbol) %>%
  mutate(Gene_sum = n()) %>% # Count peptides per gene
  ungroup() %>%
  mutate(Rank = row_number())

write_xlsx(final_scores, out_score_file)

cat("==================================================\n")
cat("Fib-SCORE Pipeline finished! Integrated scores saved to:\n", out_score_file, "\n")