# -----------------------------------------------------------------------------
# Visium HD - build CellChat objects restricted to the tumour regions
#
# Purpose   : Attach the tumour / non-tumour barcode annotation to the GP0-GP3
#             sections and rebuild the CellChat objects on tumour spots only.
#
# Inputs    : $PROJECT_ROOT/visium_hd/raw_data/GP2.rds, GP3.rds, GP0.rds
#             $PROJECT_ROOT/visium_hd/cloupe_csv/S0439254_Tumor.csv
#             $PROJECT_ROOT/visium_hd/cloupe_csv/S0441838_Tumor.csv
#             $PROJECT_ROOT/visium_hd/cloupe_csv/S0442225_Tumor.csv
#
# Outputs   : $PROJECT_ROOT/visium_hd/cellchat/GP3_tumor_cellchat.rds
#             $PROJECT_ROOT/visium_hd/cellchat/GP2_tumor_cellchat.rds
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 02_build_cellchat_tumor_regions.R
# -----------------------------------------------------------------------------

library(Seurat)
library(CellChat)
library(qs)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")

# CellChat reference database (secreted signalling)
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")

# ---- Add tumour-region barcode information (must match all cells) -----------
GP3_tumor_barcode <- read.csv(file.path(HD_ROOT, "cloupe_csv", "S0439254_Tumor.csv"))
GP3_tumor_barcodes <- GP3_tumor_barcode$Barcode[GP3_tumor_barcode$Tumor == "Tumor"]
GP3_tumor_barcodes <- sub("-1$", "", GP3_tumor_barcodes)
GP3_clean_names <- sub("-S0439254_9$", "", Cells(Seurat_GP3))
Seurat_GP3$Tumor_group <- ifelse(GP3_clean_names %in% GP3_tumor_barcodes, "Tumor", "Others")

GP2_tumor_barcode <- read.csv(file.path(HD_ROOT, "cloupe_csv", "S0441838_Tumor.csv"))
GP2_tumor_barcodes <- GP2_tumor_barcode$Barcode[GP2_tumor_barcode$Tumor == "Tumor"]
GP2_tumor_barcodes <- sub("-1$", "", GP2_tumor_barcodes)
GP2_clean_names <- sub("-S0441838_5$", "", Cells(Seurat_GP2))
Seurat_GP2$Tumor_group <- ifelse(GP2_clean_names %in% GP2_tumor_barcodes, "Tumor", "Others")

GP0_tumor_barcode <- read.csv(file.path(HD_ROOT, "cloupe_csv", "S0442225_Tumor.csv"))
GP0_tumor_barcodes <- GP0_tumor_barcode$Barcode[GP0_tumor_barcode$Tumor == "Tumor"]
GP0_tumor_barcodes <- sub("-1$", "", GP0_tumor_barcodes)
GP0_clean_names <- sub("-S0442225_10$", "", Cells(Seurat_GP0))
Seurat_GP0$Tumor_group <- ifelse(GP0_clean_names %in% GP0_tumor_barcodes, "Tumor", "Others")

# =============================================================================
# [GP3] SOLID-growth tumour region
# =============================================================================
Seurat_GP3_tumor <- subset(Seurat_GP3, subset = Tumor_group == "Tumor")
# spatial locations of the spots from the full-resolution images (first two columns only)
GP3_tumor_spatial.locs <- Seurat::GetTissueCoordinates(
  Seurat_GP3_tumor, scale = NULL, cols = c("imagerow", "imagecol"))
GP3_tumor_spatial.locs <- as.matrix(GP3_tumor_spatial.locs[, 1:2])
rownames(GP3_tumor_spatial.locs) <- Cells(Seurat_GP3_tumor)
# scale factors and spot diameters of the full-resolution images
Seurat_GP3_tumor@images[[1]]@scale.factors
GP3_tumor_scale.factors <- list(
  spot.diameter = 16, # HD 16 um
  spot = Seurat_GP3_tumor@images[[1]]@scale.factors$spot,
  fiducial = Seurat_GP3_tumor@images[[1]]@scale.factors$fiducial,
  hires = Seurat_GP3_tumor@images[[1]]@scale.factors$hires,
  lowres = Seurat_GP3_tumor@images[[1]]@scale.factors$lowres
)
spot.size <- 16 # 10X HD: spot diameter ~16 um
# conversion factor (pixel -> micrometre)
conversion.factor <- spot.size / GP3_tumor_scale.factors$spot
spatial.factors <- data.frame(
 ratio = conversion.factor,
 tol  = spot.size / 2
)
GP3_tumor_cellchat <- createCellChat(
  object = Seurat_GP3_tumor, group.by = "cell_subtype",
  assay = "SCT", datatype = "spatial",
  coordinates = GP3_tumor_spatial.locs,
  do.sparse = T,
  spatial.factors = spatial.factors)
GP3_tumor_cellchat@images$scale.factors <- GP3_tumor_scale.factors

GP3_tumor_cellchat@DB <- CellChatDB.use
GP3_tumor_cellchat <- subsetData(GP3_tumor_cellchat)
future::plan("multisession", workers = 40) # multithreading
options(future.globals.maxSize = 100000 * 1024^5)
GP3_tumor_cellchat <- identifyOverExpressedGenes(GP3_tumor_cellchat)
GP3_tumor_cellchat <- identifyOverExpressedInteractions(GP3_tumor_cellchat)
#GP3_tumor_cellchat <- projectData(GP3_tumor_cellchat, PPI.human)

GP3_tumor_cellchat <- computeCommunProb(GP3_tumor_cellchat, type="truncatedMean", trim = 0.1,
distance.use = TRUE, interaction.range = 250, scale.distance = 0.05,
contact.dependent = TRUE, contact.range = 30)
GP3_tumor_cellchat <- filterCommunication(GP3_tumor_cellchat, min.cells = 5)
GP3_tumor_cellchat <- computeCommunProbPathway(GP3_tumor_cellchat)
GP3_tumor_cellchat <- aggregateNet(GP3_tumor_cellchat)
GP3_tumor_cellchat <- netAnalysis_computeCentrality(
  GP3_tumor_cellchat, slot.name = "netP")
qsave(GP3_tumor_cellchat, file = file.path(HD_ROOT, "cellchat", "GP3_tumor_cellchat.rds"))
df.net.GP3 <- subsetCommunication(GP3_tumor_cellchat)

# =============================================================================
# [GP2] GGO-growth tumour region
# =============================================================================
Seurat_GP2_tumor <- subset(Seurat_GP2, subset = Tumor_group == "Tumor")
# spatial locations of the spots from the full-resolution images (first two columns only)
GP2_tumor_spatial.locs <- Seurat::GetTissueCoordinates(
  Seurat_GP2_tumor, scale = NULL, cols = c("imagerow", "imagecol"))
GP2_tumor_spatial.locs <- as.matrix(GP2_tumor_spatial.locs[, 1:2])
rownames(GP2_tumor_spatial.locs) <- Cells(Seurat_GP2_tumor)
# scale factors and spot diameters of the full-resolution images
Seurat_GP2_tumor@images[[1]]@scale.factors
GP2_tumor_scale.factors <- list(
  spot.diameter = 16, # HD 16 um
  spot = Seurat_GP2_tumor@images[[1]]@scale.factors$spot,
  fiducial = Seurat_GP2_tumor@images[[1]]@scale.factors$fiducial,
  hires = Seurat_GP2_tumor@images[[1]]@scale.factors$hires,
  lowres = Seurat_GP2_tumor@images[[1]]@scale.factors$lowres
)
spot.size <- 16 # 10X HD: spot diameter ~16 um
conversion.factor <- spot.size / GP2_tumor_scale.factors$spot
spatial.factors <- data.frame(
 ratio = conversion.factor,
 tol  = spot.size / 2
)
GP2_tumor_cellchat <- createCellChat(
  object = Seurat_GP2_tumor, group.by = "cell_subtype",
  assay = "SCT", datatype = "spatial",
  coordinates = GP2_tumor_spatial.locs,
  do.sparse = T,
  spatial.factors = spatial.factors)
GP2_tumor_cellchat@images$scale.factors <- GP2_tumor_scale.factors

GP2_tumor_cellchat@DB <- CellChatDB.use
GP2_tumor_cellchat <- subsetData(GP2_tumor_cellchat)
future::plan("multisession", workers = 40) # multithreading
options(future.globals.maxSize = 100000 * 1024^5)
GP2_tumor_cellchat <- identifyOverExpressedGenes(GP2_tumor_cellchat)
GP2_tumor_cellchat <- identifyOverExpressedInteractions(GP2_tumor_cellchat)
#GP2_tumor_cellchat <- projectData(GP2_tumor_cellchat, PPI.human)

GP2_tumor_cellchat <- computeCommunProb(GP2_tumor_cellchat, type="truncatedMean", trim = 0.1,
distance.use = TRUE, interaction.range = 250, scale.distance = 0.07,
contact.dependent = TRUE, contact.range = 30)
GP2_tumor_cellchat <- filterCommunication(GP2_tumor_cellchat, min.cells = 10) # ~25% removed
GP2_tumor_cellchat <- computeCommunProbPathway(GP2_tumor_cellchat)
GP2_tumor_cellchat <- aggregateNet(GP2_tumor_cellchat)
GP2_tumor_cellchat <- netAnalysis_computeCentrality(
  GP2_tumor_cellchat, slot.name = "netP", future.seed = TRUE)
saveRDS(GP2_tumor_cellchat, file = file.path(HD_ROOT, "cellchat", "GP2_tumor_cellchat.rds"))
df.net.GP2 <- subsetCommunication(GP2_tumor_cellchat)
