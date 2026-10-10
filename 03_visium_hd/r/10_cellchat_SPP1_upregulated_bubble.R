# -----------------------------------------------------------------------------
# Visium HD - SPP1+ macrophage outgoing signalling bubble
#
# Purpose   : Draw the CellChat bubble plot of the FN1 and SPP1 signalling
#             pathways sent by SPP1+ macrophages (ma_SPP1) to the other cell
#             subtypes, comparing the two tumour-region datasets.
#
# Inputs    : cellchat.rds (CellChat merged object; produced upstream)
#
# Outputs   : netVisual_bubble_SPP1.png and netVisual_bubble_SPP1.pdf
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 10_cellchat_SPP1_upregulated_bubble.R
# -----------------------------------------------------------------------------

library(CellChat)
library(patchwork)
library(glue)
library(ComplexHeatmap)
library(optparse)
library(circlize)

cellchat <- readRDS("cellchat.rds")

p1 <- netVisual_bubble(cellchat, sources.use = "ma_SPP1", targets.use = c("apCAF", "AT2", "B_naive", "ma_MARCO", "ma_SPP1",
                "malignant_CCL28", "malignant_CXCL14", "malignant_myoCAF",
                "malignant_plasma", "plasma", "SDP", "T_eff", "T_ex"), comparison = c(1, 2), angle.x = 45, line.on = FALSE,
                 signaling = c("FN1", "SPP1"))
p1

ggsave(p1, file = "netVisual_bubble_SPP1.png", dpi = 300, width = 15, height = 5)
ggsave(p1, file = "netVisual_bubble_SPP1.pdf", width = 15, height = 5)
