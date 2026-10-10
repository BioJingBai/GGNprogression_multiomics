# -----------------------------------------------------------------------------
# Setup, cell-type annotation and myeloid/lymphoid sub-clustering (CosMx SMI)
#
# Purpose   : Load the CosMx SMI Seurat object and the sample sheet, define the
#             study grouping variables (Group7 / group / tissue), build the
#             final Level-2 / Level-3 / Level-4 cell-type labels, and
#             sub-cluster macrophages, T cells and B cells so that the
#             SPP1+ / CCL18+ / STAB1+ / MARCO+ / COX2+ macrophage states and the
#             CD4T / CD8T / NKT / B / B_GC states are available downstream.
#
# Inputs    : <PROJECT_ROOT>/raw_data/seurat.rds        (CosMx SMI Seurat object)
#             <PROJECT_ROOT>/raw_data/sample_info.csv   (sample -> Group7 map)
#
# Outputs   : <PROJECT_ROOT>/processed/smi_annotated.qs            (annotated object)
#             <PROJECT_ROOT>/outputs/Macrophages_cluster_markers_SCT.csv
#             <PROJECT_ROOT>/outputs/Tcell_cluster_markers_SCT_0.3.csv
#             <PROJECT_ROOT>/outputs/Bcell_cluster_markers_SCT_0.5.csv
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 00_setup_and_annotation.R
# -----------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Dependencies
## ---------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(stringr)
library(qs)
library(conflicted)

# Set the project root; no absolute machine paths are hard-coded.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
dir.create(file.path(PROJECT_ROOT, "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PROJECT_ROOT, "outputs"),   recursive = TRUE, showWarnings = FALSE)


## ---------------------------------------------------------------------------
## 1. Load data
## ---------------------------------------------------------------------------
# Read the SMI Seurat object and the sample sheet
smi <- readRDS(file.path(PROJECT_ROOT, "raw_data", "seurat.rds"))
sample_info <- read.csv(file.path(PROJECT_ROOT, "raw_data", "sample_info.csv"))
sample_info[1, 3] <- "1-9"
sample_info[2, 3] <- "10-18"


## ---------------------------------------------------------------------------
## 2. Global options
## ---------------------------------------------------------------------------
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::arrange)
conflicted::conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(base::union)
set.seed(1234)


## ---------------------------------------------------------------------------
## 3. Add the study grouping (Group7 / group / tissue) to the metadata
## ---------------------------------------------------------------------------
# NOTE: this block operates on `smi1`, which is not defined yet; the object
# loaded above is `smi`. Uncomment the next line to run the block from scratch.
# smi1 <- smi

# Make the Sample column and smi@meta.data$sample the same type
sample_info$Sample <- as.character(sample_info$Sample)
sample_info$Group7 <- as.character(sample_info$Group7)
smi1@meta.data$sample <- as.character(smi1@meta.data$sample)
# Create the Group7 column by matching sample to Group7
smi1@meta.data$Group7 <- sample_info$Group7[match(smi1@meta.data$sample, sample_info$Sample)]
smi1@meta.data$Group7 <- as.factor(smi1@meta.data$Group7)
# Check the result
head(smi1@meta.data[, c("sample", "Group7")])
# Cells that do not belong to GP5 (drop the GP5_T / GP5_TE / GP5_N samples)
cells_to_keep <- rownames(smi1@meta.data[!(smi1@meta.data$Group7 %in% c("GP5_T", "GP5_TE", "GP5_N")), ])
# Subset using the cells argument
smi1 <- subset(smi1, cells = cells_to_keep)
# Groups that remain
unique(smi1@meta.data$Group7)
# Derive group and tissue from Group7
smi1@meta.data$group <- sub("_.*", "", smi1@meta.data$Group7)
smi1@meta.data$tissue <- sub(".*_", "", smi1@meta.data$Group7)
smi1@meta.data$group <- as.factor(smi1@meta.data$group)
smi1@meta.data$tissue <- as.factor(smi1@meta.data$tissue)


## ---------------------------------------------------------------------------
## 4. Final Level-2 cell type
## ---------------------------------------------------------------------------
unique(smi@meta.data$Level2_celltype)
unique(smi@meta.data$Level3_celltype)

smi$Level2_celltype <- as.character(smi$Level2_celltype)
smi$Level3_celltype <- as.character(smi$Level3_celltype)
smi$final_Level2_celltype <- NA
smi$final_Level2_celltype[smi$Level2_celltype %in% c(
  "T", "B", "Plasma", "Macrophages", "DC",
  "NoFs", "CAFs", "mesothelial", "Endothelial", "SMC", "EP", "Mix"
)] <- smi$Level2_celltype[smi$Level2_celltype %in% c(
  "T", "B", "Plasma", "Macrophages", "DC",
  "NoFs", "CAFs", "mesothelial", "Endothelial", "SMC", "EP", "Mix"
)]
smi$final_Level2_celltype[smi$Level3_celltype %in% c(
  "CD4T", "CD8T", "NKT"
)] <- "T"
smi$final_Level2_celltype[smi$Level3_celltype %in% c(
  "activated_B", "Naive_B"
)] <- "B"
smi$final_Level2_celltype[smi$Level3_celltype %in% c(
  "Epithelial_MUC1", "Epithelial_PTEN", "Epithelial_F2R",
  "Epithelial_IFI30", "Epithelial_CHMP2B"
)] <- "EP"
smi$final_Level2_celltype[smi$Level3_celltype %in% c(
  "endothelia_Bro", "endothelia_Lym"
)] <- "Endo"
smi$final_Level2_celltype[smi$Level3_celltype %in% c(
  "Mix"
)] <- "Mix"
smi$final_Level2_celltype[smi$Level3_celltype %in% c(
  "Plasma"
)] <- "Plasma"
smi$final_Level2_celltype[smi$Level3_celltype %in% c("CAFs")] <- "CAF"
smi$final_Level2_celltype[smi$Level3_celltype %in% c("NoFs")] <- "NoF"
smi$final_Level2_celltype[smi$Level3_celltype %in% c("mesothelial")] <- "mesothelial"
table(smi$final_Level2_celltype, useNA = "ifany")

smi$final_Level2_celltype <- factor(smi$final_Level2_celltype,
    levels = c("EP", "T", "B", "Plasma", "Macrophages", "DC", "Endo", "CAF", "NoF", "mesothelial", "SMC", "Mix")
)


## ---------------------------------------------------------------------------
## 5. Re-cluster macrophages (ma_COX2 / ma_MARCO / ma_CCL18 / ma_STAB1 / ma_SPP1)
## ---------------------------------------------------------------------------
macrophages <- subset(smi, subset = final_Level2_celltype == "Macrophages")
macrophages <- RunPCA(macrophages, assay = "SCT", features = VariableFeatures(macrophages))
# Build the neighbour graph on the SCT assay
macrophages <- FindNeighbors(macrophages, dims = 1:20, graph.name = "SCT_snn", assay = "SCT")
macrophages <- FindClusters(macrophages, resolution = 0.5, graph.name = "SCT_snn")
macrophages <- RunUMAP(macrophages, dims = 1:20)
DimPlot(macrophages, reduction = "umap", label = TRUE)
macrophage_markers <- FindAllMarkers(macrophages,
                                     assay = "SCT",
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25)
write.csv(macrophage_markers,
          file = file.path(PROJECT_ROOT, "outputs", "Macrophages_cluster_markers_SCT.csv"),
          row.names = FALSE)

cluster_names <- c(
  "0" = "ma_COX2",   # COX1 COX2 SLC11A1
  "1" = "ma_COX2",
  "2" = "ma_MARCO",
  "3" = "ma_CCL18",
  "4" = "ma_STAB1",
  "5" = "ma_SPP1",
  "6" = "ma_MARCO",
  "7" = "Plasma"
)
# Build a vector as long as the number of cells
cluster_annotation <- cluster_names[as.character(macrophages$seurat_clusters)]
names(cluster_annotation) <- Cells(macrophages)
macrophages <- AddMetaData(macrophages, metadata = cluster_annotation, col.name = "Macrophage_cluster_name")


## ---------------------------------------------------------------------------
## 6. Re-cluster T cells (CD4T / CD8T / NKT / T_cm / T_TIGHT / T_eff)
## ---------------------------------------------------------------------------
Tcell <- subset(smi, subset = final_Level2_celltype == "T")
Tcell <- RunPCA(Tcell, assay = "SCT", features = VariableFeatures(Tcell))
# Build the neighbour graph on the SCT assay
Tcell <- FindNeighbors(Tcell, dims = 1:20, graph.name = "SCT_snn", assay = "SCT")
Tcell <- FindClusters(Tcell, resolution = 0.5, graph.name = "SCT_snn")
Tcell <- RunUMAP(Tcell, dims = 1:20)
DimPlot(Tcell, reduction = "umap", label = TRUE)
Tcell_markers <- FindAllMarkers(Tcell,
                                     assay = "SCT",
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25)
write.csv(Tcell_markers,
          file = file.path(PROJECT_ROOT, "outputs", "Tcell_cluster_markers_SCT_0.3.csv"),
          row.names = FALSE)

cluster_names <- c(
  "0" = "EP",
  "1" = "T",
  "2" = "T_TIGHT",
  "3" = "T_cm",     # IL7R
  "4" = "T",        # LTB CD74
  "5" = "B",        # BCL2
  "6" = "T_eff",    # CXCR4
  "7" = "T",
  "8" = "T",
  "9" = "Plasma"
)
cluster_annotation <- cluster_names[as.character(Tcell$seurat_clusters)]
names(cluster_annotation) <- Cells(Tcell)
Tcell <- AddMetaData(Tcell, metadata = cluster_annotation, col.name = "Tcell_cluster_name")


## ---------------------------------------------------------------------------
## 7. Re-cluster B cells (B / B_GC)
## ---------------------------------------------------------------------------
Bcell <- subset(smi, subset = final_Level2_celltype == "B")
Bcell <- RunPCA(Bcell, assay = "SCT", features = VariableFeatures(Bcell))
# Build the neighbour graph on the SCT assay
Bcell <- FindNeighbors(Bcell, dims = 1:20, graph.name = "SCT_snn", assay = "SCT")
Bcell <- FindClusters(Bcell, resolution = 0.5, graph.name = "SCT_snn")
Bcell <- RunUMAP(Bcell, dims = 1:20)
DimPlot(Bcell, reduction = "umap", label = TRUE)
Bcell_markers <- FindAllMarkers(Bcell,
                                     assay = "SCT",
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25)
write.csv(Bcell_markers,
          file = file.path(PROJECT_ROOT, "outputs", "Bcell_cluster_markers_SCT_0.5.csv"),
          row.names = FALSE)

genes_to_check <- c("CD21", "CD23")
genes_to_check <- str_to_upper(genes_to_check)  # convert to upper case

p <- DotPlot(Bcell, features = unique(genes_to_check), assay = "RNA") + coord_flip()
p

cluster_names <- c(
  "0" = "B",
  "1" = "B_GC",     # CXCL13 CXCR4
  "2" = "B",
  "3" = "B",        # LTB
  "4" = "B_GC",     # CXCL13 CXCR4
  "5" = "B",        # BLK
  "6" = "B",        # CD37
  "7" = "B",        # PAX5
  "8" = "B",
  "9" = "Plasma"
)
cluster_annotation <- cluster_names[as.character(Bcell$seurat_clusters)]
names(cluster_annotation) <- Cells(Bcell)
Bcell <- AddMetaData(Bcell, metadata = cluster_annotation, col.name = "Bcell_cluster_name")


## ---------------------------------------------------------------------------
## 8. Merge the T / B / macrophage sub-clusters back into final_Level3_celltype
## ---------------------------------------------------------------------------
# Create the new columns in the large Seurat object first
smi$final_Level2_celltype <- as.character(smi$final_Level2_celltype)
Tcell$Tcell_cluster_name <- as.character(Tcell$Tcell_cluster_name)
Bcell$Bcell_cluster_name <- as.character(Bcell$Bcell_cluster_name)
macrophages$Macrophage_cluster_name <- as.character(macrophages$Macrophage_cluster_name)
smi$final_Level3_celltype <- smi$final_Level2_celltype
# Use match to keep the indices aligned
smi$final_Level3_celltype[match(Cells(Tcell), Cells(smi))] <- Tcell$Tcell_cluster_name
smi$final_Level3_celltype[match(Cells(Bcell), Cells(smi))] <- Bcell$Bcell_cluster_name
smi$final_Level3_celltype[match(Cells(macrophages), Cells(smi))] <- macrophages$Macrophage_cluster_name
# Check the result
table(smi$final_Level3_celltype, useNA = "ifany")
smi$final_Level2_celltype <- factor(smi$final_Level2_celltype,
    levels = c("EP", "T", "B", "Plasma", "Macrophages", "DC", "Endo", "CAF", "NoF", "mesothelial", "SMC", "Mix")
)
smi$final_Level3_celltype <- factor(smi$final_Level3_celltype,
    levels = c("EP", "T", "T_cm", "T_TIGHT", "B", "B_GC", "Plasma", "ma_SPP1", "ma_CCL18", "ma_STAB1", "ma_HLA-DR", "ma_MARCO", "ma_COX2", "DC", "Endo", "CAF", "NoF", "mesothelial", "SMC", "Mix")
)


## ---------------------------------------------------------------------------
## 9. Final Level-4 cell type (re-cluster only the macrophages)
## ---------------------------------------------------------------------------
# Keep the CD4 / CD8 / activated / naive labels of Level3 and only refresh the
# macrophage layer with the sub-clustering above.
smi$Level3_celltype <- as.character(smi$Level3_celltype)
macrophages$Macrophage_cluster_name <- as.character(macrophages$Macrophage_cluster_name)
smi$final_Level4_celltype <- smi$Level3_celltype
smi$final_Level4_celltype[match(Cells(macrophages), Cells(smi))] <- macrophages$Macrophage_cluster_name
smi$final_Level4_celltype <- gsub("^Epithelial.*", "EP", smi$final_Level4_celltype)
# Check the result
table(smi$final_Level4_celltype, useNA = "ifany")
smi$final_Level2_celltype <- factor(smi$final_Level2_celltype,
    levels = c("EP", "T", "B", "Plasma", "Macrophages", "DC", "Endo", "CAF", "NoF", "mesothelial", "SMC", "Mix")
)

# Re-create the column after the environment was cleared in the source session
smi$Level3_celltype <- as.character(smi$Level3_celltype)
smi$final_Level3_celltype <- as.character(smi$final_Level3_celltype)
smi$final_Level4_celltype <- smi$Level3_celltype
macrophage_cells <- smi$Level3_celltype == "Macrophages"
smi$final_Level4_celltype[macrophage_cells] <- smi$final_Level3_celltype[macrophage_cells]
smi$final_Level4_celltype <- gsub("^Epithelial_.*", "EP", smi$final_Level4_celltype)
smi$final_Level4_celltype[smi$final_Level4_celltype == "ma_HLA-DR"] <- "ma_MARCO"
smi$final_Level4_celltype <- factor(smi$final_Level4_celltype,
    levels = c("EP", "CD4T", "CD8T", "NKT", "activated_B", "Naive_B", "Plasma", "ma_SPP1", "ma_CCL18", "ma_STAB1", "ma_MARCO", "ma_COX2", "DC", "endothelia_Bro", "endothelia_Lym", "CAFs", "NoFs", "mesothelial", "SMC", "Mix"))
table(smi$final_Level4_celltype, useNA = "ifany")


## ---------------------------------------------------------------------------
## 10. Save the annotated object for the downstream scripts
## ---------------------------------------------------------------------------
qsave(smi, file.path(PROJECT_ROOT, "processed", "smi_annotated.qs"))
