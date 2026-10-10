# -----------------------------------------------------------------------------
# Visium HD - CCL19-CCR7 interactions among lymphocyte subtypes (GGO-growth)
#
# Purpose   : Show the CCL communication between the naive-B / effector-T /
#             Tfh compartments in the GGO-growth section (GP2 / S0441838_5),
#             as a CellChat ligand-receptor bubble plot, and add the CCL19 /
#             CCR7 expression violins per cell subtype for the same section.
#
# Inputs    : $PROJECT_ROOT/visium_hd/external_cellchat/single_sample_min50/tumor/S0441838_5/cellchat.rds
#             $PROJECT_ROOT/visium_hd/raw_data/subtype_reannotation.rData
#
# Outputs   : $PROJECT_ROOT/visium_hd/figures/CCL_bubble.pdf
#             $PROJECT_ROOT/visium_hd/figures/CCL19_GP2.pdf
#             $PROJECT_ROOT/visium_hd/figures/CCR7_GP2.pdf
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 07_cellchat_CCL19_CCR7_lymphocytes.R
# -----------------------------------------------------------------------------

library(Seurat)
library(CellChat)
library(ggplot2)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")
dir.create(file.path(HD_ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)

# ---- CCL19-CCR7 bubble among lymphocyte subtypes (GGO-growth) ---------------
P_cellchat <- readRDS(file.path(HD_ROOT, "external_cellchat",
                                "single_sample_min50", "tumor",
                                "S0441838_5", "cellchat.rds"))
target_pathway <- "CCL"
# unique(P_cellchat@netP$pathways) or unique(P_cellchat@net$pathways)
from_cells <- c("B_naive", "T_eff", "Tfh")
to_cells <- c("T_eff", "Tfh")
df <- subsetCommunication(P_cellchat, signaling = target_pathway)
df_filtered <- df %>%
  dplyr::filter(source %in% from_cells & target %in% to_cells)
p <- netVisual_bubble(P_cellchat, signaling = target_pathway, sources.use = from_cells, targets.use = to_cells)
ggsave(file.path(HD_ROOT, "figures", "CCL_bubble.pdf"), p, width = 4, height = 2)

# ---- CCL19 and CCR7 expression per cell subtype (GP2) -----------------------
load(file.path(HD_ROOT, "raw_data", "subtype_reannotation.rData"))
# cell types kept for the violin plots
target_subtypes <- c(
  # T cell subtypes
  "T_ex", "T_eff", "T_proApoptosis", "Tfh", "T_naive", "CD4T_APC",
  # B cell subtypes
  "B_mem", "B_active", "plasma", "B_naive",
  # epithelial subtypes
  "malignant_CXCL14", "malignant_CCL28", "AT1", "AT2", "club", "ciliated",
  # macrophage subtypes
  "ma_MARCO", "ma_NR4A1", "ma_FOLR2", "ma_EPAS1", "ma_SPP1", "cDC",
  # CAF subtypes
  "proCAF", "myoCAF", "apCAF", "Fib_ADH1B", "iCAF"
)
selected_cells <- object@meta.data[object@meta.data$cell_subtype %in% target_subtypes, ]
object_filtered <- subset(object, cells = rownames(selected_cells))
object <- subset(object_filtered, Group %in% c("S0441838_5"))
true_cells <- which(object@assays$Spatial@cells@`.Data`[, "counts"])
object@assays$Spatial@cells <- object@assays$Spatial@cells[true_cells, ]
cluster_colors <- c("#EA7070", "#FDC4B6", "#E59572", "#2694AB", "#96CEB4", "#FFEEAD", "#D9534F",
                    "#FFAD60", "#57D1C9", "#ED5485","#7F95D1", "#FFE869", "#DE4307", "#F29C2B",
                    "#F6D04D", "#8BC24C", "#4695D6", "#FED95C", "#FA6E57", "#F69E53",
                    "#248888", "#E6E6E6","#E7475E", "#F0D879", "#29A2C6", "#1f77b4", "#ff7f0e")

GP2 <- object
p2 <- VlnPlot(GP2, features = "CCL19", group.by = "cell_subtype", cols = cluster_colors, pt.size = 0) +
  ggtitle("CCL19_GP2") +
  theme(axis.title.x = element_blank())  # drop the x-axis title
p2
ggsave(p2, file = file.path(HD_ROOT, "figures", "CCL19_GP2.pdf"), width = 12, height = 12)

p2 <- VlnPlot(GP2, features = "CCR7", group.by = "cell_subtype", cols = cluster_colors, pt.size = 0) +
  ggtitle("CCR7_GP2") +
  theme(axis.title.x = element_blank())
p2
ggsave(p2, file = file.path(HD_ROOT, "figures", "CCR7_GP2.pdf"), width = 12, height = 12)
