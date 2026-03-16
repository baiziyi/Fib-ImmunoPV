# ==============================================================================
# scRNA-seq Heatmap, K-means Clustering & GO Enrichment Pipeline
# ==============================================================================

rm(list=ls()) # Clean environment

# Load libraries
library(tidyverse)
library(Seurat)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)   

# ==============================================================================
# 1. Global Parameters
# ==============================================================================

# --- 1.1 Input Settings ---
sc_obj            <- scRNA_dataset             # Seurat object
target_features   <- target_geneset$Symbol     # Target gene set
group_name        <- "Ctr"                     # Output prefix
celltype_col      <- "celltype"                # Metadata cell type column
species_org_db    <- org.Mm.eg.db              # org.db

# --- 1.2 Clustering Settings ---
k_clusters        <- 5

# --- 1.3 Target Cell Types (Ordered for heatmap) ---
target_celltypes  <- c(
  "b cell", "plasma cell", "CD4T", "CD8T", "nk cell",
  "neutrophil", "basophil", "plasmacytoid dendritic cell", "dendritic cell",
  "classical monocyte","intermediate monocyte","non-classical monocyte", "macrophage",
  "smooth muscle cell", "bronchial smooth muscle cell", "vascular associated smooth muscle cell",
  "pericyte cell", "mesothelial cell", "adventitial cell", "fibroblast", "alveolar fibroblast", "myofibroblast cell",
  "type i pneumocyte", "type ii pneumocyte", "respiratory mucous cell", "respiratory goblet cell",
  "lung ciliated cell", "club cell", "basal cell", "pulmonary ionocyte", "serous cell of epithelium of bronchus",
  "endothelial cell of artery", "vein endothelial cell", "capillary endothelial cell", 
  "capillary aerocyte", "lung microvascular endothelial cell", 
  "Bronchial vessel endothelial cell", "endothelial cell of lymphatic vessel"
)

# Cell type color dictionary
color_dict <- c(
  'b cell'='#E64B35FF', 'plasma cell'='#efefef', 'CD4T'= '#002240', 
  'CD8T'= 'blue', 'nk cell'='#7aa6dc', 'neutrophil'='#00A087FF',
  'basophil'='#3C5488FF', 'plasmacytoid dendritic cell'='#F39B7FFF',
  'dendritic cell'='#8491B4FF', 'classical monocyte'='#91D1C2FF',
  'intermediate monocyte'="#0099B4FF", 'non-classical monocyte'="#d69794",
  'macrophage'="red", 'smooth muscle cell'="#42B540FF",
  'bronchial smooth muscle cell'='#eed78c', 'vascular associated smooth muscle cell'='#cc7564',
  'pericyte cell'='#80c291', 'mesothelial cell'='#79a2e8', 'adventitial cell'="#c76ea7",
  'fibroblast'="#d4b598", 'alveolar fibroblast'="#e3edc3", 'myofibroblast cell'='#f38f8c',
  'type i pneumocyte'='#ecb871', 'type ii pneumocyte'='#eebed6',
  'respiratory mucous cell'='#7E6148FF', 'respiratory goblet cell'="#F39B7FFF",
  'lung ciliated cell'="#8fc4ef", 'club cell'='#E64B35FF', 'basal cell'="#d69794",
  'pulmonary ionocyte'="#e1e588", 'serous cell of epithelium of bronchus'="#a1c6d3",
  'endothelial cell of artery'="#e2aca6", 'vein endothelial cell'="#d3faae",
  'capillary endothelial cell'="#abb993", 'capillary aerocyte'="#f28cb6",
  'lung microvascular endothelial cell'="#cfaef4", 'Bronchial vessel endothelial cell'="#95d1a1",
  'endothelial cell of lymphatic vessel'="#f7f998"
)

# --- 1.4 Output Paths ---
out_dir_cluster   <- "/home/Project/Lung_fibrosis_Immupep/Figure_Data/Kmeans_cluster"
out_dir_enrich    <- "/home/Project/Lung_fibrosis_Immupep/Figure_Data/Kmeans_enrichment"

# Create dirs safely
dir.create(out_dir_cluster, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir_enrich, group_name), recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. Matrix Extraction
# ==============================================================================
cat("\n[1/4] Calculating average expression...\n")

# Get average expression
avg_exp_obj <- AverageExpression(
  sc_obj, 
  assays = "RNA", 
  group.by = celltype_col, 
  features = target_features, 
  return.seurat = FALSE
)$RNA

# Filter available cell types
available_cells <- intersect(target_celltypes, colnames(avg_exp_obj))
exp_matrix      <- as.matrix(avg_exp_obj[, available_cells])


# ==============================================================================
# 3. Heatmap & K-means
# ==============================================================================
cat("\n[2/4] Generating Pheatmap...\n")

# Column annotations
annotation_col <- data.frame(Group = rep(group_name, length(available_cells)), Type = available_cells)
rownames(annotation_col) <- available_cells

# Annotation colors
anno_colors <- list(Type = color_dict[available_cells])

# Generate heatmap (silent)
pep_heatmap <- pheatmap(
  exp_matrix, 
  kmeans_k          = k_clusters, 
  border            = FALSE, 
  annotation_col    = annotation_col, 
  annotation_colors = anno_colors,
  cluster_cols      = TRUE,
  silent            = TRUE 
)

# Save outputs
heatmap_pdf_path <- file.path(out_dir_cluster, paste0(group_name, "_cluster.pdf"))
pdf(heatmap_pdf_path, width = 14, height = 5)
print(pep_heatmap)
dev.off()

saveRDS(pep_heatmap, file.path(out_dir_cluster, paste0(group_name, "_cluster.rds")))
cat("      Heatmap cluster saved to:", heatmap_pdf_path, "\n")


# ==============================================================================
# 4. GO Enrichment per Cluster
# ==============================================================================
cat("\n[3/4] Running GO Enrichment...\n")

cluster_df <- data.frame(cluster = pep_heatmap$kmeans$cluster)
all_enrich_results <- list()

for(i in 1:k_clusters) {
  cat(sprintf("  -> Processing Cluster %d...\n", i))
  
  # Get cluster genes
  genes_in_cluster <- rownames(cluster_df)[cluster_df$cluster == i]
  
  # Run enrichGO
  ego <- enrichGO(
    gene          = genes_in_cluster,
    OrgDb         = species_org_db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  # Simplify and save valid results
  if(!is.null(ego) && nrow(ego@result) > 0) {
    ego_sim <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    
    res_df <- as.data.frame(ego_sim) %>% mutate(cluster = as.character(i))
    all_enrich_results[[as.character(i)]] <- res_df
    
    out_file <- file.path(out_dir_enrich, group_name, paste0("Cluster_", i, "_GO_0.7.xlsx"))
    writexl::write_xlsx(res_df, out_file)
  }
}


# ==============================================================================
# 5. Visualization (Top GO Terms)
# ==============================================================================
cat("\n[4/4] Generating Dot Plot...\n")

# Combine results
combined_enrich_df <- bind_rows(all_enrich_results)

# Select Top 5 GO terms per cluster
plot_df <- combined_enrich_df %>%
  group_by(cluster) %>%
  slice_min(order_by = qvalue, n = 5) %>% 
  ungroup()

if(nrow(plot_df) > 0) {
  # Format GeneRatio
  plot_df$GeneRatio <- sapply(plot_df$GeneRatio, function(x) {
    nums <- as.numeric(unlist(strsplit(x, "/")))
    nums[1] / nums[2]
  })
  
  max_log_q <- max(-log10(plot_df$qvalue))
  
  # Plot tile + point
  p_enrich <- ggplot(plot_df, aes(x = cluster, y = str_wrap(Description, width = 45), fill = -log10(qvalue))) +
    geom_tile() +
    geom_point(aes(size = Count), color = "orange") +
    scale_fill_gradient(low = "white", high = "#6C80A2", guide = "colorbar", limits = c(0, max_log_q)) + 
    labs(x = "K-means Cluster", y = NULL, title = paste(group_name, "Enrichment")) +
    theme_bw() +
    theme(
      panel.grid   = element_blank(),
      panel.border = element_blank(),
      axis.ticks   = element_blank(),
      axis.text.y  = element_text(size = 10, color = "black"),
      axis.text.x  = element_text(size = 12, color = "black", face = "bold"),
      plot.title   = element_text(size = 15, hjust = 0.5)
    ) +
    scale_y_discrete(position = "right")
  
  # Save plot
  out_pdf <- file.path(out_dir_enrich, paste0(group_name, "_Enrich_Dotplot.pdf"))
  pdf(file = out_pdf, width = 6.8, height = 5.5)
  print(p_enrich)
  dev.off()
  
  # Export merged results
  writexl::write_xlsx(combined_enrich_df, file.path(out_dir_enrich, paste0(group_name, "_All_Enrichment.xlsx")))
  
  cat("      Plot saved to:", out_pdf, "\n")
  cat("==================================================\n")
} else {
  cat("[WARNING] No significant GO terms found.\n")
}