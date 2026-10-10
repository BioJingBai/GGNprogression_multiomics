# -----------------------------------------------------------------------------
# Visium HD - compare tumour-region communication, SOLID-growth vs GGO-growth
#
# Purpose   : Merge the two tumour-region CellChat objects (GP2 = GGO-growth,
#             GP3 = SOLID-growth) and compare their interaction counts /
#             strengths, differential interaction networks, signalling-role
#             scatter plots and functional network similarity, to identify the
#             cell-cell interactions that are up- or down-regulated in
#             SOLID-growth relative to GGO-growth.
#
# Inputs    : $PROJECT_ROOT/visium_hd/cellchat/GP2_tumor_cellchat.rds
#             $PROJECT_ROOT/visium_hd/cellchat/GP3_tumor_cellchat.rds
#
# Outputs   : $PROJECT_ROOT/visium_hd/figures/compare_GP23_tumor*.pdf
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 12_cellchat_compare_GGO_vs_SOLID_tumor.R
# -----------------------------------------------------------------------------

library(CellChat)
library(Seurat)
library(patchwork)
library(qs)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")
fig_dir      <- file.path(HD_ROOT, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

GP2_tumor_cellchat <- readRDS(file.path(HD_ROOT, "cellchat", "GP2_tumor_cellchat.rds"))
GP3_tumor_cellchat <- readRDS(file.path(HD_ROOT, "cellchat", "GP3_tumor_cellchat.rds"))

setwd(fig_dir)
object.list <- list(GP2 = GP2_tumor_cellchat, GP3 = GP3_tumor_cellchat)
tumor_cellchat <- mergeCellChat(object.list, add.names = names(object.list))
pdf("compare_GP23_tumor.pdf", height = 8, width = 15)
gg1 <- compareInteractions(tumor_cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(tumor_cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

pdf("compare_GP23_tumor_cirecle.pdf", height = 8, width = 15)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(tumor_cellchat, weight.scale = T)
netVisual_diffInteraction(tumor_cellchat, weight.scale = T, measure = "weight")
dev.off()

pdf("compare_GP23_tumor_heatmap.pdf", height = 8, width = 15)
gg1 <- netVisual_heatmap(tumor_cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(tumor_cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

# With multiple datasets one can display the interaction count / strength between
# any two cell groups per dataset. To control node size and edge weight across
# datasets, compute the maximum cell number per group and the maximum interaction
# count (or weight) over all datasets.
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# cell groups with significant signalling changes between datasets
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
pdf("compare_GP23_tumor_signal_celltype.pdf", height = 8, width = 15)
patchwork::wrap_plots(plots = gg)
dev.off()

# compare structural or functional differences
tumor_cellchat <- computeNetSimilarityPairwise(tumor_cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
tumor_cellchat <- netEmbedding(tumor_cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
tumor_cellchat <- netClustering(tumor_cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualisation in 2D space
netVisual_embeddingPairwise(tumor_cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
