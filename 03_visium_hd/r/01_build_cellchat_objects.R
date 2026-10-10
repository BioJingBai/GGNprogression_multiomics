# -----------------------------------------------------------------------------
# Visium HD - build the whole-section CellChat objects (GP3 and GP2)
#
# Purpose   : Standardise each section with SCTransform and infer spatially
#             constrained cell-cell communication with CellChat from the spot
#             coordinates.
#
# Inputs    : $PROJECT_ROOT/visium_hd/raw_data/GP2.rds, GP3.rds, GP0.rds, GP1.rds
#
# Outputs   : $PROJECT_ROOT/visium_hd/cellchat/GP3_cellchat.rds
#             $PROJECT_ROOT/visium_hd/cellchat/GP2_cellchat.rds
#             (SCTransform-normalised GP0-GP3.rds written back to raw_data/)
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 01_build_cellchat_objects.R
# -----------------------------------------------------------------------------

library(Seurat)
library(CellChat)
library(qs)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")

# ---- CellChat reference database -------------------------------------------
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse
# subset the database to the secreted-signalling fraction
CellChatDB.use <- subsetDB(CellChatDB,
                           search = "Secreted Signaling", # Cell-Cell Contact, ECM-Receptor, Non-protein Signaling, Secreted Signaling
                           key = "annotation")

# ---- Normalisation (SCTransform + log-normalise) ----------------------------
# the two normalisations gave the same result; the latter is used downstream
Seurat_GP3 <- SCTransform(Seurat_GP3, assay = "Spatial", verbose = FALSE)
Seurat_GP2 <- SCTransform(Seurat_GP2, assay = "Spatial", verbose = FALSE)
Seurat_GP1 <- SCTransform(Seurat_GP1, assay = "Spatial", verbose = FALSE)
Seurat_GP0 <- SCTransform(Seurat_GP0, assay = "Spatial", verbose = FALSE)
qsave(Seurat_GP0, file.path(HD_ROOT, "raw_data", "GP0.rds"))
qsave(Seurat_GP1, file.path(HD_ROOT, "raw_data", "GP1.rds"))
qsave(Seurat_GP2, file.path(HD_ROOT, "raw_data", "GP2.rds"))
qsave(Seurat_GP3, file.path(HD_ROOT, "raw_data", "GP3.rds"))

# =============================================================================
# [GP3] SOLID-growth section
# =============================================================================
# spatial locations of the spots from the full-resolution images (first two columns only)
GP3_spatial.locs <- Seurat::GetTissueCoordinates(
  Seurat_GP3, scale = NULL, cols = c("imagerow", "imagecol"))
GP3_spatial.locs <- as.matrix(GP3_spatial.locs[, 1:2])
rownames(GP3_spatial.locs) <- Cells(Seurat_GP3)
# scale factors and spot diameters of the full-resolution images
Seurat_GP3@images[[1]]@scale.factors
GP3_scale.factors <- list(
  spot.diameter = 16, # HD 16 um
  spot = Seurat_GP3@images[[1]]@scale.factors$spot,
  fiducial = Seurat_GP3@images[[1]]@scale.factors$fiducial,
  hires = Seurat_GP3@images[[1]]@scale.factors$hires,
  lowres = Seurat_GP3@images[[1]]@scale.factors$lowres
)
spot.size <- 16 # 10X HD: spot diameter ~16 um
# conversion factor (pixel -> micrometre)
conversion.factor <- spot.size / GP3_scale.factors$spot # the denominator is the length in pixels
spatial.factors <- data.frame(
 ratio = conversion.factor,
 tol  = spot.size / 2 # tolerance = half a spot
)
# Create a CellChat object. Note: scale.factors has to be attached separately,
# it cannot be passed directly to the constructor.
GP3_cellchat <- createCellChat(
  object = Seurat_GP3, group.by = "cell_subtype",
  assay = "SCT", datatype = "spatial",
  coordinates = GP3_spatial.locs,
  do.sparse = T,
  spatial.factors = spatial.factors)
GP3_cellchat@images$scale.factors <- GP3_scale.factors

# Set the ligand-receptor interaction database
GP3_cellchat@DB <- CellChatDB.use
# Preprocess the expression data for communication analysis: subset the
# expression data of the signalling genes
GP3_cellchat <- subsetData(GP3_cellchat)
future::plan("multisession", workers = 40) # multithreading
options(future.globals.maxSize = 100000 * 1024^5)
# identify over-expressed genes
GP3_cellchat <- identifyOverExpressedGenes(GP3_cellchat)
# identify over-expressed ligand-receptor pairs
GP3_cellchat <- identifyOverExpressedInteractions(GP3_cellchat)
# project the gene expression data onto the PPI
#GP3_cellchat <- projectData(GP3_cellchat, PPI.human)

# Compute the communication probability and infer the cellular communication
# network. (The earlier function does not work here; the um-based version is
# used because the pixel-to-um conversion has already been applied.)
#GP3_cellchat <- computeCommunProb(GP3_cellchat,
  #type = "truncatedMean", # mean gene expression per cell group
  #trim = 0.1,
  #distance.use = TRUE,   # spatial distance constraint on communication probability
  #scale.distance =0.07)
GP3_cellchat <- computeCommunProb(GP3_cellchat, type="truncatedMean", trim = 0.07,
distance.use = TRUE, interaction.range = 250, scale.distance = 0.05,
contact.dependent = TRUE, contact.range = 30)
# filter out cell-cell communication when a group contains only few cells
GP3_cellchat <- filterCommunication(GP3_cellchat, min.cells = 10) # ~25% removed
# infer communication at the signalling-pathway level
GP3_cellchat <- computeCommunProbPathway(GP3_cellchat)
# compute the aggregated cell-cell communication network
GP3_cellchat <- aggregateNet(GP3_cellchat)

# extract the information of all communication pairs
df.net <- subsetCommunication(GP3_cellchat)

# inspect the available pathways
unique(df.net$pathway_name)

# keep the SPP1 pathway only
df.SPP1 <- df.net[df.net$pathway_name == "SPP1", ]

# inspect SPP1-CD44 communication
df.SPP1[df.SPP1$receptor == "CD44", ]
# compute the network centrality score
GP3_cellchat <- netAnalysis_computeCentrality(
  GP3_cellchat, slot.name = "netP") # "netP" = inferred cell-cell communication network of signalling pathways
# save the CellChat object
saveRDS(GP3_cellchat, file = file.path(HD_ROOT, "cellchat", "GP3_cellchat.rds"))

# =============================================================================
# [GP2] GGO-growth section
# =============================================================================
# spatial locations of the spots from the full-resolution images (first two columns only)
GP2_spatial.locs <- Seurat::GetTissueCoordinates(
  Seurat_GP2, scale = NULL, cols = c("imagerow", "imagecol"))
GP2_spatial.locs <- as.matrix(GP2_spatial.locs[, 1:2])
rownames(GP2_spatial.locs) <- Cells(Seurat_GP2)
# scale factors and spot diameters of the full-resolution images
Seurat_GP2@images[[1]]@scale.factors
GP2_scale.factors <- list(
  spot.diameter = 16, # HD 16 um
  spot = Seurat_GP2@images[[1]]@scale.factors$spot,
  fiducial = Seurat_GP2@images[[1]]@scale.factors$fiducial,
  hires = Seurat_GP2@images[[1]]@scale.factors$hires,
  lowres = Seurat_GP2@images[[1]]@scale.factors$lowres
)
spot.size <- 16 # 10X HD: spot diameter ~16 um
# conversion factor (pixel -> micrometre)
conversion.factor <- spot.size / GP2_scale.factors$spot # the denominator is the length in pixels
spatial.factors <- data.frame(
 ratio = conversion.factor,
 tol  = spot.size / 2 # tolerance = half a spot
)
# Create a CellChat object (scale.factors attached separately)
GP2_cellchat <- createCellChat(
  object = Seurat_GP2, group.by = "cell_subtype",
  assay = "SCT", datatype = "spatial",
  coordinates = GP2_spatial.locs,
  do.sparse = T,
  spatial.factors = spatial.factors)
GP2_cellchat@images$scale.factors <- GP2_scale.factors

# Set the ligand-receptor interaction database
GP2_cellchat@DB <- CellChatDB.use
# preprocess the expression data for communication analysis
GP2_cellchat <- subsetData(GP2_cellchat)
future::plan("multisession", workers = 40) # multithreading
options(future.globals.maxSize = 100000 * 1024^5)
GP2_cellchat <- identifyOverExpressedGenes(GP2_cellchat)
GP2_cellchat <- identifyOverExpressedInteractions(GP2_cellchat)
#GP2_cellchat <- projectData(GP2_cellchat, PPI.human)

GP2_cellchat <- computeCommunProb(GP2_cellchat, type="truncatedMean", trim = 0.1,
distance.use = TRUE, interaction.range = 250, scale.distance = 0.05,
contact.dependent = TRUE, contact.range = 30)
GP2_cellchat <- filterCommunication(GP2_cellchat, min.cells = 10) # ~25% removed
GP2_cellchat <- computeCommunProbPathway(GP2_cellchat)
GP2_cellchat <- aggregateNet(GP2_cellchat)
GP2_cellchat <- netAnalysis_computeCentrality(
  GP2_cellchat, slot.name = "netP")
saveRDS(GP2_cellchat, file = file.path(HD_ROOT, "cellchat", "GP2_cellchat.rds"))
