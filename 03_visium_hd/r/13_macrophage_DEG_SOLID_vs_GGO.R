# -----------------------------------------------------------------------------
# Visium HD - differentially expressed genes in macrophages, SOLID vs GGO
#
# Purpose   : Find the genes that are up- or down-regulated in macrophages of
#             the SOLID-growth section (S0439254_9) relative to the GGO-growth
#             section (S0441838_5).
#
# Inputs    : $PROJECT_ROOT/visium_hd/raw_data/All_Seurat.rds
#
# Outputs   : DE_genes / up_genes_GP3 / down_genes_GP3 data frames in the
#             workspace.
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 13_macrophage_DEG_SOLID_vs_GGO.R
# -----------------------------------------------------------------------------

library(Seurat)
library(qs)
library(dplyr)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")

All_Seurat <- qread(file.path(HD_ROOT, "raw_data", "All_Seurat.rds"))

# Assume Level2_celltype is in meta.data and Sample is in meta.data
macrophages <- subset(All_Seurat, subset = cell_type_renew == "Macrophage")
# Wilcoxon test by default
DE_genes <- FindMarkers(
  object = macrophages,
  ident.1 = "S0439254_9",
  ident.2 = "S0441838_5",
  group.by = "Group",   # grouping column
  logfc.threshold = 0.25, # adjustable threshold
  min.pct = 0.1           # expressed in at least 10% of cells
)

up_genes_GP3 <- DE_genes %>%
  dplyr::filter(avg_log2FC > 0 & p_val_adj < 0.05) %>%
  dplyr::arrange(desc(avg_log2FC))


down_genes_GP3 <- DE_genes %>%
  dplyr::filter(avg_log2FC < 0 & p_val_adj < 0.05) %>%
  dplyr::arrange(desc(avg_log2FC))
