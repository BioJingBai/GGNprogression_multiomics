# -----------------------------------------------------------------------------
# Visium HD - environment setup, niche annotation and four-section split
#
# Purpose   : Load the packages, colour palettes and Seurat objects used by every
#             downstream Visium HD script; attach the BANKSY niche clustering to
#             the merged object; harmonise the cell-subtype labels; and split the
#             merged object into the four sequenced sections (GP0-GP3).
#
# Inputs    : $PROJECT_ROOT/visium_hd/raw_data/subtype_reannotation.rData      (annotated Seurat object)
#             $PROJECT_ROOT/visium_hd/raw_data/integrate_BANKSY_Plot_simp.rData (BANKSY result)
#             $PROJECT_ROOT/visium_hd/raw_data/All_Seurat.rds
#             $PROJECT_ROOT/visium_hd/raw_data/GP0.rds ... GP3.rds
#             $PROJECT_ROOT/visium_hd/cellchat/GP2_cellchat.rds       (pre-built CellChat object)
#
# Outputs   : refreshed GP0-GP3.rds / All_Seurat.rds (qsave) under raw_data/
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 00_setup_and_annotation.R
# -----------------------------------------------------------------------------

# ---- Load packages ----------------------------------------------------------
options(stringsAsFactors = F)
library(RColorBrewer)
library(wesanderson)
library(paletteer)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)
library(ggpubr)
library(plyr)
library(ggbeeswarm)
library(argparser)
#library(scales)
library(SingleCellExperiment)
library(SpatialExperiment)
library(imcRtools)
library(diffcyt)
library(cowplot)
library(dittoSeq)
library(viridis)
library(bluster)
library(BiocParallel)
library(tidyverse)
library(scater)
library(harmony)
library(FastPG)
library(glue)

library(readxl)
library(DemographicTable)
library(tidyverse)
library(writexl)
library(conflicted)
library(flextable)
library(officer)
library(dplyr)
library(stats)
library(ggplot2)
library(ggpubr)
library(Biobase)
library(data.table)
library(org.Mm.eg.db)
library(clustree)

library(ggrepel)
library(corrplot)
library(cowplot)
library(ComplexHeatmap)
library(vcd)
library(psych)
library(monocle)
library(tidyr)
library(tibble)
library(readr)
library(janitor)
library(reshape)
library(survminer)
library(broom)
library(ggstatsplot)
library(hrbrthemes)
library(ggnewscale)
library(stringr)

library(MatrixGenerics)
library(circlize)
library(gridBase)
library(clusterProfiler)
library(GOplot)
library(openxlsx)
library(msigdbr)
library(GSEABase)
library(GSVA)
library(pheatmap)
library(RColorBrewer)
library(enrichplot)
library(EnrichedHeatmap)
library(limma)

library(survival)
library(glmnet)
library(ranger)
library(pROC)
library(survivalROC)

library(Seurat)
library(hdf5r)
library(data.table)
library(patchwork)
#if (!require("remotes")) install.packages("remotes")
#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)
library(harmony)
library(RColorBrewer)
library(EnhancedVolcano)
library(COSG)
library(scRNAtoolVis)
library(CellChat)
library(qs)
library(monocle)
library(biomaRt)
library(homologene)

# install a large-language-model single-cell annotation helper
#devtools::install_github("cafferychen777/mLLMCelltype", subdir="R")
library(mLLMCelltype)
library(cowplot) # used for plot_grid
# SCpubr for publication-grade visualisation
#if (!requireNamespace("SCpubr", quietly = TRUE)) {remotes::install_github("enblacar/SCpubr")}
library(SCpubr)

# ---- Colour configuration ---------------------------------------------------
pdf("colorRampPalette_brewer_pal.pdf", width = 8, height = 10)
display.brewer.all() # inspect available palettes
dev.off()

mycol <- c("#E64B35FF","#4DBBD5FF",
        "#F5E4BB",
        "#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF", "#DC0000FF","#7E6148FF","#B09C85FF","gray",brewer.pal(12,"Paired"),brewer.pal(6,"Accent")[6],"#666666","#B2182B",brewer.pal(8,"Set2"))
MYCOLOR <- c("#FF6700", "#9370DB", "#F8D568","#00AD43", "#89CFF0", "#BA160C",           "#FF91AF", "#A6A6A6", "#006DB0","#C154C1", "#D99A6C", "#96C8A2",           "#FBEC5D", "#77B5FE", "#E29CD2","#007F5C", "#ACBF60", "#7B68EE",           "#00FFCD", "#D3AF37", "#50C878","#1974D2", "#FF0000", "#9F00FF",           "#91A3B0", "#8B4513", "#4166F5","#C19A6B", "#6EAEA1", "#39FF14",           "#0247FE", "#AF6E4D", "#FF66CC","#9400D3", "#A57164", "#00FFFF",           "#DF00FF", "#FFAA1D", "#20B2AA","#CB6D51", "#5218FA", "#BC8F8F","#000000")
MYCOLOR2 <- c(  "#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF",  "#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF",  "#d1cc78","#b2387e","#E64B354C","#4DBBD54C",  "#00A0874C","#3C54884C","#F39B7F4C","#8491B44C",  "#91D1C24C","#DC00004C","#fcec7c","#ff92e0",  "#E64B357F","#4DBBD57F","#00A0877F","#3C54887F",  "#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F",  "#eae8c1","#ef4db1","#E64B35B2","#4DBBD5B2",  "#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2",  "#91D1C2B2","#DC0000B2","#edd600","#c64ea4")

# ---- Load data --------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")

load(file.path(HD_ROOT, "raw_data", "subtype_reannotation.rData"))
All_Seurat <- object
load(file.path(HD_ROOT, "raw_data", "integrate_BANKSY_Plot_simp.rData"))
banksy <- object
remove(object)
All_Seurat <- qread(file.path(HD_ROOT, "raw_data", "All_Seurat.rds"))
Seurat_GP3 <- qread(file.path(HD_ROOT, "raw_data", "GP3.rds"))
Seurat_GP2 <- qread(file.path(HD_ROOT, "raw_data", "GP2.rds"))
Seurat_GP1 <- qread(file.path(HD_ROOT, "raw_data", "GP1.rds"))
Seurat_GP0 <- qread(file.path(HD_ROOT, "raw_data", "GP0.rds"))
GP2_cellchat <- readRDS(file.path(HD_ROOT, "cellchat", "GP2_cellchat.rds"))

# ---- Save data --------------------------------------------------------------
qsave(Seurat_GP0, file.path(HD_ROOT, "raw_data", "GP0.rds"))
qsave(Seurat_GP2, file.path(HD_ROOT, "raw_data", "GP2.rds"))
qsave(Seurat_GP3, file.path(HD_ROOT, "raw_data", "GP3.rds"))
qsave(Seurat_GP1, file.path(HD_ROOT, "raw_data", "GP1.rds"))
qsave(All_Seurat, file.path(HD_ROOT, "raw_data", "All_Seurat.rds"))

# ---- Niche analysis ---------------------------------------------------------
# 1. add the BANKSY niche clustering to the metadata
All_Seurat <- AddMetaData(All_Seurat, metadata = banksy@meta.data)
All_Seurat@graphs$BANKSY_nn <- banksy@graphs$BANKSY_nn
All_Seurat@graphs$BANKSY_snn <- banksy@graphs$BANKSY_snn
All_Seurat@reductions$pca.banksy <- banksy@reductions$pca.banksy
All_Seurat@reductions$harmony.banksy <- banksy@reductions$harmony
All_Seurat@reductions$umap.banksy <- banksy@reductions$umap
# rename the two labels that were carried over from the IMC nomenclature
All_Seurat@meta.data$cell_subtype <- All_Seurat@meta.data$cell_subtype %>%
  dplyr::recode(
    "endo" = "Endothelial",
    "SMC_MYH11" = "SMCs"
  )
# extract the four sequenced sections
Seurat_GP0 <- subset(All_Seurat, subset = orig.ident == "S0442225_10")
Seurat_GP1 <- subset(All_Seurat, subset = orig.ident == "S0432519_9")
Seurat_GP2 <- subset(All_Seurat, subset = orig.ident == "S0441838_5")
Seurat_GP3 <- subset(All_Seurat, subset = orig.ident == "S0439254_9")
