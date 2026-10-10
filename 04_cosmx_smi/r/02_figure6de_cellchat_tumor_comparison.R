# -----------------------------------------------------------------------------
# CellChat comparison of solid (MP) versus part-solid (P) tumours
#
# Purpose   : Merge the P_T and MP_T spatial CellChat objects and compare the
#             two groups: overall interaction number/strength, differential
#             cell-cell signalling, SPP1 pathway bubble plots, and the
#             incoming/outgoing pathway-strength heatmaps.
#
# Inputs    : <PROJECT_ROOT>/processed/cellchat/smi_P_T_cellchat202511_updatecell2.0.qs
#             <PROJECT_ROOT>/processed/cellchat/smi_MP_T_cellchat202511_updatecell2.0.qs
#
# Outputs   : <PROJECT_ROOT>/processed/cellchat/smi_cellchat_Tumor202511_updatecell2.0.qs
#             <PROJECT_ROOT>/figures/cellchat_compare/Tumor_cellchat_*.pdf
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 02_figure6de_cellchat_tumor_comparison.R
# -----------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Dependencies
## ---------------------------------------------------------------------------
library(Seurat)
library(CellChat)
library(ComplexHeatmap)
library(patchwork)
library(ggplot2)
library(qs)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
dir.create(file.path(PROJECT_ROOT, "figures", "cellchat_compare"), recursive = TRUE, showWarnings = FALSE)

## Load the CellChat objects produced by 01_cellchat_build_final_level4.R
smi_P_T_cellchat  <- qread(file.path(PROJECT_ROOT, "processed", "cellchat", "smi_P_T_cellchat202511_updatecell2.0.qs"))
smi_MP_T_cellchat <- qread(file.path(PROJECT_ROOT, "processed", "cellchat", "smi_MP_T_cellchat202511_updatecell2.0.qs"))


## ---------------------------------------------------------------------------
## 3. Tumour comparison figures
## ---------------------------------------------------------------------------
# SPP1 also highlighted CXCL12-CXCR4 as a potential route, as it is enriched in MP
object.list <- list(P = smi_P_T_cellchat, MP = smi_MP_T_cellchat)
cellchat_Tumor <- mergeCellChat(object.list, add.names = names(object.list))
qsave(cellchat_Tumor, file.path(PROJECT_ROOT, "processed", "cellchat", "smi_cellchat_Tumor202511_updatecell2.0.qs"))
# Overall interaction comparison
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_completenumber.pdf"), width = 15, height = 8)
gg1 <- compareInteractions(cellchat_Tumor, show.legend = F,
                           group = c(1, 2))
gg2 <- compareInteractions(cellchat_Tumor, show.legend = F,
                           group = c(1, 2), measure = "weight")
gg1 + gg2
dev.off()

# Number and strength of interactions between cell types — circle plot
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_complete_number.pdf"), width = 15, height = 8)
netVisual_diffInteraction(cellchat_Tumor, weight.scale = T, measure = "weight")
dev.off()
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_complete_strength.pdf"), width = 15, height = 8)
netVisual_diffInteraction(cellchat_Tumor, weight.scale = T)
dev.off()

# Number and strength of interactions between cell types — heatmap
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_complete_heatmap.pdf"), width = 15, height = 8)
gg1 <- netVisual_heatmap(cellchat_Tumor)
gg2 <- netVisual_heatmap(cellchat_Tumor, measure = "weight")
gg1 + gg2
dev.off()

# Number and strength of interactions between selected cell types — circle plot
# Grouped according to celltype_order
# `celltype_order` is not used here; see the header note
# Grouping used only for a compact visual summary; it carries no biological meaning
group.cellType <- c(rep("T/NK", 9), rep("Group1", 4), rep("Group2", 5))
group.cellType <- factor(group.cellType, levels = c("T/NK", "Group1", "Group2"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat_Tumor <- mergeCellChat(object.list, add.names = names(object.list))
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents", "count", "count.merged"))
# interactions
par(mfrow = c(1, 2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge = T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
# weight.merged
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat_Tumor, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat_Tumor, weight.scale = T, measure = "weight.merged", label.edge = T)


# Cells with a significant difference in the signals they send or receive
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))  # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
ggsave(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffcell.pdf"), width = 25, height = 12)


# Identify the signalling changes of specific cell groups
# Visualise how the outgoing and incoming signals change from P to MP
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffsignal.pdf"), width = 22, height = 10)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, comparison = c(1, 2), idents.use = "ma_SPP1", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "ma_STAB1", signaling.exclude = "MIF")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "ma_CCL18", signaling.exclude = "MIF")
gg4 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "ma_MARCO", signaling.exclude = "MIF")
gg5 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "ma_COX2", signaling.exclude = "MIF")
gg6 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "CD4T", signaling.exclude = "MIF")
gg7 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "CD8T", signaling.exclude = "MIF")
gg8 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "NKT", signaling.exclude = "MIF")  # MIF really dominates the difference
gg9 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "activated_B", signaling.exclude = "MIF")
gg10 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "Naive_B", signaling.exclude = "MIF")
gg11 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "Plasma", signaling.exclude = "MIF")
gg12 <- netAnalysis_signalingChanges_scatter(cellchat_Tumor, idents.use = "DC", signaling.exclude = "MIF")
patchwork::wrap_plots(plots = list(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, gg10, gg11))
dev.off()

## Identify the different signalling pathways
# Compare the overall information flow of each pathway or ligand-receptor pair
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffpathway.pdf"), width = 15, height = 8)
gg1 <- rankNet(cellchat_Tumor, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_Tumor, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()

# Compare the outgoing and incoming signalling patterns of each cell group
i = 1
# Combine all identified signalling pathways from the different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i + 1], width = 8, height = 16)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 8, height = 16)
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i + 1], width = 8, height = 16)
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffcellsignal_outgoing.pdf"), width = 10, height = 18)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffcellsignal_incoming.pdf"), width = 10, height = 18)
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))
dev.off()


## Identify the up- and down-regulated ligand-receptor pairs
levels(cellchat_Tumor@idents$joint)
# Identify dysfunctional signalling by comparing communication probabilities
gg1 <- netVisual_bubble(cellchat_Tumor, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in MP", angle.x = 45, remove.isolate = T, signaling = "CCL")
gg2 <- netVisual_bubble(cellchat_Tumor, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MP", angle.x = 45, remove.isolate = T, signaling = "CCL")
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffcellLR_CCL.pdf"), width = 20, height = 14)
gg1 + gg2
dev.off()

gg3 <- netVisual_bubble(cellchat_Tumor, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in MP", angle.x = 45, remove.isolate = T, signaling = "CXCL")
gg4 <- netVisual_bubble(cellchat_Tumor, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MP", angle.x = 45, remove.isolate = T, signaling = "CXCL")
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffcellLR_CXCL.pdf"), width = 20, height = 14)
gg3 + gg4
dev.off()

gg5 <- netVisual_bubble(cellchat_Tumor, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in MP", angle.x = 45, remove.isolate = T, signaling = "SPP1")
gg6 <- netVisual_bubble(cellchat_Tumor, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MP", angle.x = 45, remove.isolate = T, signaling = "SPP1")
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffcellLR_SPP1.pdf"), width = 20, height = 14)
gg5 + gg6
dev.off()

gg7 <- netVisual_bubble(cellchat_Tumor, sources.use = c(2:4, 11, 12), targets.use = c(9, 10), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in MP", angle.x = 45, remove.isolate = T)
gg8 <- netVisual_bubble(cellchat_Tumor, sources.use = c(2:4, 11, 12), targets.use = c(9, 10), comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MP", angle.x = 45, remove.isolate = T)
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffcellLR_maccaf_B.pdf"), width = 12, height = 12)
gg7 + gg8
dev.off()

gg9 <- netVisual_bubble(cellchat_Tumor, sources.use = c(2:4, 11, 12), targets.use = c(6, 7, 8), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in MP", angle.x = 45, remove.isolate = T)
gg10 <- netVisual_bubble(cellchat_Tumor, sources.use = c(2:4, 11, 12), targets.use = c(6, 7, 8), comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MP", angle.x = 45, remove.isolate = T)
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffcellLR_maccaf_T.pdf"), width = 15, height = 12)
gg9 + gg10
dev.off()

gg11 <- netVisual_bubble(cellchat_Tumor, sources.use = c(9, 10), targets.use = c(6, 7, 8), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in MP", angle.x = 45, remove.isolate = T)
gg12 <- netVisual_bubble(cellchat_Tumor, sources.use = c(9, 10), targets.use = c(6, 7, 8), comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MP", angle.x = 45, remove.isolate = T)
pdf(file.path(PROJECT_ROOT, "figures", "cellchat_compare", "Tumor_cellchat_diffcellLR_B_T.pdf"), width = 12, height = 12)
gg11 + gg12
dev.off()

# Identify dysfunctional signalling by differential expression analysis
# (the block above used communication probabilities instead)
# Set the "experimental group"
pos.dataset = "P"
features.name = paste0(pos.dataset, ".merged")
cellchat_Tumor <- identifyOverExpressedGenes(cellchat_Tumor, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05, thresh.p = 0.05, group.DE.combined = FALSE)
net <- netMappingDEG(cellchat_Tumor, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat_Tumor, net = net, datasets = "left", ligand.logFC = NULL, receptor.logFC = 0.05)
net.down <- subsetCommunication(cellchat_Tumor, net = net, datasets = "right", ligand.logFC = NULL, receptor.logFC = -0.05)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_Tumor)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_Tumor)
# Find all significant outgoing / incoming / both signals for the chosen features and cell groups
df <- findEnrichedSignaling(object.list[[2]], features = c("CCL19", "SPP1"), pattern = "outgoing")

## Visualise the identified up- / down-regulated ligand-receptor pairs
# Bubble plot
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_Tumor, pairLR.use = pairLR.use.up, sources.use = 16, targets.use = c(1:15), comparison = c(1, 2), angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat_Tumor, pairLR.use = pairLR.use.down, sources.use = 16, targets.use = c(1:15), comparison = c(1, 2), angle.x = 90, remove.isolate = T, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# Chord diagram
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = "net", net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = "net", net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# Word cloud
library(wordcloud)
computeEnrichmentScore(net.down, species = "human", variable.both = TRUE)
computeEnrichmentScore(net.up, species = "human", variable.both = TRUE)

## Compare the receptors of a pathway
cellchat_Tumor@meta$datasets = factor(cellchat_Tumor@meta$datasets, levels = c("P", "MP"))  # set factor level
plotGeneExpression(cellchat_Tumor, signaling = "SPP1", split.by = "datasets", colors.ggplot = T)
