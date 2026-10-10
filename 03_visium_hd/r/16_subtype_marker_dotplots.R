# -----------------------------------------------------------------------------
# Visium HD - marker dotplots for every cell subtype
#
# Purpose   : Draw one jjDotPlot per compartment (T, B, CAF, epithelial,
#             myeloid) showing the curated marker genes across the annotated
#             subtypes, using the spatially-resolved "Spatial" assay.
#
# Inputs    : $PROJECT_ROOT/visium_hd/subannots/{T_cell,B_cell,CAFs,Epi,Macrophage}/all_Anno.rData
#             $PROJECT_ROOT/visium_hd/subannots/T_cell/marker_T.txt
#
# Outputs   : Marker_DotPlot_{T,B,CAF,Epi,macro}.pdf in the working directory
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 16_subtype_marker_dotplots.R
# -----------------------------------------------------------------------------

# load packages
library(Seurat)
#library(tidydr)
#library(scplotter)
#library(scCustomize)
#library(SCP)
library(tidyverse)
library(ggpubr)
library(ggalluvial)
library(cowplot)
library(scRNAtoolVis)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")
annot_dir    <- file.path(HD_ROOT, "subannots")

# ==================== T cell subtypes ====================
load(file.path(annot_dir, "T_cell", "all_Anno.rData"))
all_genes <- c("TNFRSF4", "LAG3", "TIGIT", "GZMA", "ISG15", "DKK3",
               "CD3E", "CD3D", "CD8A", "TCF7", "CCR7", "LEF1", "SELL", "GZMK",
               "CXCL13", "PDCD1", "NR4A1", "KLRG1", "IL7R",
               "BAG3", "CASP10", "BAK1", "TXNIP", "FAS", "KLRG1",
               "CD4", "TCF7", "SELL", "NKG7", "CXCL13", "CXCR4",
               "CD8A", "CCR7", "CCR10", "CTLA4",
               "CD4", "TCF7", "KLRG1")

# remove duplicated genes
all_genes_unique <- unique(all_genes)
all_genes_unique <- as.data.frame(all_genes_unique)
all_genes_unique
head(object)
# build the gene data frame
markerGene <- data.frame(gene = all_genes_unique)
markerGene$gene <- markerGene$all_genes_unique

P1 <- jjDotPlot(object = object, gene = markerGene$gene,
          assay = "Spatial",
          anno = T, id = "cell_sub_type",
          #plot.margin = c(5, 3, 3, 1),
          ytree = F,
          same.pos.label = T,
          textRot = 0,
          rescale = T,
          x.text.angle = 90,
          rescale.min = -1,
          rescale.max = 3,
          dot.col = c('blue','yellow','red'))
P1
ggsave("Marker_DotPlot_T.pdf", P1, width = 10, height = 4)
table(object@meta.data$cell_sub_type)
jjDotPlot(object = object,
          gene = markerGene$gene,
          anno = T,
          plot.margin = c(3,1,1,1),
          tree.pos = 'left',
          same.pos.label = T,
          yPosition = 10.3)
table(object@meta.data$SCT_snn_res.0.5)

# read the marker list to visualise
Marker_table <- read.table(file.path(annot_dir, "T_cell", "marker_T.txt"), sep = "\t", header = T)
cell.order <- colnames(Marker_table)
#Marker_table %>% as.list %>% unlist %>% unname()
Marker <- Marker_table %>%
  pivot_longer(cols = everything(), names_to = "cluster", values_to = "gene") %>% filter(gene != "")
Marker <- Marker[Marker$gene %in% rownames(object),]
Marker$cluster <- factor(Marker$cluster, levels = level_order)
Marker <- Marker[!duplicated(Marker$gene),]
# subset to the T cell subtypes directly
object@meta.data$Celltype <- object@meta.data$cell_sub_type

object <- subset(object,
                 subset = Celltype %in% c("T_ex", "T_eff", "T_proApoptosis", "Tfh", "T_naive", "CD4T_APC"))

# set the display order of the cell types
level_order <- c("T_ex", "T_eff", "T_proApoptosis", "Tfh", "T_naive", "CD4T_APC")
object@meta.data$Celltype <- factor(object@meta.data$Celltype, levels = level_order)

# set the cell identities
Idents(object) <- "Celltype"

# check the result
print("cell type distribution after subsetting:")
print(table(object$Celltype))
print("display order:")
print(levels(object$Celltype))
print(table(object$SCT_snn_res.0.5))

#object <- subset(ref, BANKSY_snn_res.0.2 %in% c("0", "4", "6"))
true_cells <- which(object@assays$Spatial@cells@`.Data`[,"counts"])
object@assays$Spatial@cells <- object@assays$Spatial@cells[true_cells,]
levels(Marker$cluster)
P1 <- jjDotPlot(object = object, markerGene = Marker,
          plot.margin = c(3,1,1,0),
          ytree = F,
          id = 'Celltype',
          same.pos.label = T,
          #x.text.angle = 45,
          textRot = 0,
          dot.col = c('blue','yellow','red'),
          rescale = T,
          rescale.min = -1,
          rescale.max = 3)
P1

# ==================== B cell subtypes ====================
rm(list = ls())
load(file.path(annot_dir, "B_cell", "all_Anno.rData"))
all_genes <- c("AGER", "CAV1", "S100A10", "S100A4", "LMO2", "ISG15", "CRIP2", "ITGB1",
           "EGR1", "NR4A2", "MT1X", "CD69", "CD83", "TNF",
           "XBP1", "JCHAIN", "IGHA1", "IGKC", "MZB1", "IGHG3", "NR4A2", "CD79A",
           "IGHG1", "IGLC1", "IGHM",
           "MS4A1", "CD79A", "CD79B", "CCR7", "TCL1A", "FCER2", "IGHD", "IFIT3", "ISG15")

# remove duplicated genes
all_genes_unique <- unique(all_genes)
all_genes_unique <- as.data.frame(all_genes_unique)
all_genes_unique
head(object)
table(object@meta.data$cell_sub_type)
object <- subset(object, cell_sub_type %in% c("B_mem", "B_active", "plasma", "B_naive"))
true_cells <- which(object@assays$Spatial@cells@`.Data`[,"counts"])
object@assays$Spatial@cells <- object@assays$Spatial@cells[true_cells,]
# build the gene data frame
markerGene <- data.frame(gene = all_genes_unique)
markerGene$gene <- markerGene$all_genes_unique

P1 <- jjDotPlot(object = object, gene = markerGene$gene,
          assay = "Spatial",
          anno = T, id = "cell_sub_type",
          #plot.margin = c(5, 3, 3, 1),
          ytree = F,
          same.pos.label = T,
          textRot = 0,
          rescale = T,
          x.text.angle = 90,
          rescale.min = -1,
          rescale.max = 3,
          dot.col = c('blue','yellow','red'))
P1

ggsave("Marker_DotPlot_B.pdf", P1, width = 10, height = 4)

# ==================== CAF subtypes ====================
rm(list = ls())
load(file.path(annot_dir, "CAFs", "all_Anno.rData"))
all_genes <- c("IGF1", "OGN", "C7", "COL15A1", "PI16",
           "COL10A1", "POSTN", "MMP11", "CTHRC1", "SFRP2", "LRRC15", "HOPX", "MMP1", "SFRP2",
           "ADH1B", "A2M", "C7",
           "CD74", "CXCL2", "CXCL12", "CXCL14",
           "IL6", "HGF", "HSPA6")

# remove duplicated genes
all_genes_unique <- unique(all_genes)
all_genes_unique <- as.data.frame(all_genes_unique)
all_genes_unique
head(object)
table(object@meta.data$cell_sub_type)
object <- subset(object, cell_sub_type %in% c("proCAF", "myoCAF", "Fib_ADH1B", "apCAF", "iCAF"))
true_cells <- which(object@assays$Spatial@cells@`.Data`[,"counts"])
object@assays$Spatial@cells <- object@assays$Spatial@cells[true_cells,]
# build the gene data frame
markerGene <- data.frame(gene = all_genes_unique)
markerGene$gene <- markerGene$all_genes_unique

P1 <- jjDotPlot(object = object, gene = markerGene$gene,
          assay = "Spatial",
          anno = T, id = "cell_sub_type",
          #plot.margin = c(5, 3, 3, 1),
          ytree = F,
          same.pos.label = T,
          textRot = 0,
          rescale = T,
          x.text.angle = 90,
          rescale.min = -1,
          rescale.max = 3,
          dot.col = c('blue','yellow','red'))
P1

ggsave("Marker_DotPlot_CAF.pdf", P1, width = 10, height = 4)

# ==================== epithelial subtypes ====================
rm(list = ls())
load(file.path(annot_dir, "Epi", "all_Anno.rData"))
all_genes <- c("EPCAM", "CXCL14",
           "SFTPA1", "AGER", "SFTPC", "SFTPA2", "SFTPB", "SFTPC", "WIF1",
           "EPCAM", "CCL28",
           "AGER", "CLDN18",
           "SCGB3A2", "SCGB3A1", "SCGB1A1", "SCGB1A1",
           "FOXJ1", "TPPP3", "CIMAP3")

# remove duplicated genes
all_genes_unique <- unique(all_genes)
all_genes_unique <- as.data.frame(all_genes_unique)
all_genes_unique
head(object)
table(object@meta.data$cell_sub_type)
object <- subset(object, cell_sub_type %in% c("malignant_CXCL14", "AT2", "malignant_CCL28", "AT1", "club", "ciliated"))
true_cells <- which(object@assays$Spatial@cells@`.Data`[,"counts"])
object@assays$Spatial@cells <- object@assays$Spatial@cells[true_cells,]
# build the gene data frame
markerGene <- data.frame(gene = all_genes_unique)
markerGene$gene <- markerGene$all_genes_unique

P1 <- jjDotPlot(object = object, gene = markerGene$gene,
          assay = "Spatial",
          anno = T, id = "cell_sub_type",
          #plot.margin = c(5, 3, 3, 1),
          ytree = F,
          same.pos.label = T,
          textRot = 0,
          rescale = T,
          x.text.angle = 90,
          rescale.min = -1,
          rescale.max = 3,
          dot.col = c('blue','yellow','red'))
P1

ggsave("Marker_DotPlot_Epi.pdf", P1, width = 10, height = 4)

# ==================== myeloid subtypes ====================
rm(list = ls())
load(file.path(annot_dir, "Macrophage", "all_Anno.rData"))
all_genes <- c("MARCO", "C1QB", "C1QC", "CD68",
           "IRF7", "S100A8", "NR4A1",
           "LAMP3",
           "FOLR2", "CD68",
           "LYZ", "EPAS1",
           "SPP1", "MMP9", "CSF1R")

# remove duplicated genes
all_genes_unique <- unique(all_genes)
all_genes_unique <- as.data.frame(all_genes_unique)
all_genes_unique
head(object)
table(object@meta.data$cell_sub_type)
object <- subset(object, cell_sub_type %in% c("ma_MARCO", "ma_NR4A1", "cDC", "ma_FOLR2", "ma_EPAS1", "ma_SPP1"))
true_cells <- which(object@assays$Spatial@cells@`.Data`[,"counts"])
object@assays$Spatial@cells <- object@assays$Spatial@cells[true_cells,]
# build the gene data frame
markerGene <- data.frame(gene = all_genes_unique)
markerGene$gene <- markerGene$all_genes_unique

P1 <- jjDotPlot(object = object, gene = markerGene$gene,
          assay = "Spatial",
          anno = T, id = "cell_sub_type",
          #plot.margin = c(5, 3, 3, 1),
          ytree = F,
          same.pos.label = T,
          textRot = 0,
          rescale = T,
          x.text.angle = 90,
          rescale.min = -1,
          rescale.max = 3,
          dot.col = c('blue','yellow','red'))
P1

ggsave("Marker_DotPlot_macro.pdf", P1, width = 10, height = 4)
