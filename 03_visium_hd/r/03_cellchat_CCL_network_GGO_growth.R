# -----------------------------------------------------------------------------
# Visium HD - CCL signalling network in the GGO-growth section
#
# Purpose   : Draw the CellChat CCL signalling network inferred for the
#             GGO-growth section (GP2 / S0441838_5) and the CCL ligand-receptor
#             bubble plot between the lymphocyte compartments.
#
# Inputs    : $PROJECT_ROOT/visium_hd/cellchat/GP2_cellchat.rds
#
# Outputs   : $PROJECT_ROOT/visium_hd/figures/CCL_circle_GP2.pdf (on-screen)
#             $PROJECT_ROOT/visium_hd/figures/P_T_cellchat_bubbleplot_nont.pdf
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 03_cellchat_CCL_network_GGO_growth.R
# -----------------------------------------------------------------------------

library(CellChat)
library(Seurat)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")
dir.create(file.path(HD_ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)

GP2_cellchat <- readRDS(file.path(HD_ROOT, "cellchat", "GP2_cellchat.rds"))

# ---- CCL pathway network (circle layout) ------------------------------------
netVisual_aggregate(GP2_cellchat, signaling = "CCL", layout = "circle")

# ---- CCL ligand-receptor bubble plot between lymphocyte compartments --------
levels(GP2_cellchat@idents)
netVisual_bubble(GP2_cellchat, sources.use = c(1:6, 8:12, 14:17),
                 targets.use = c(2:6), remove.isolate = FALSE)
ggsave(file.path(HD_ROOT, "figures", "P_T_cellchat_bubbleplot_nont.pdf"),
       width = 20, height = 15)
