# -----------------------------------------------------------------------------
# Spatial CellChat inference on the final Level-4 cell-type annotation
#
# Purpose   : Build the four spatial CellChat objects used for the tumour
#             comparison (P_T, P_TE, MP_T, MP_TE) on `final_Level4_celltype`,
#             restricting the ligand-receptor database to secreted signalling.
#
# Inputs    : <PROJECT_ROOT>/processed/smi_annotated.qs
#
# Outputs   : <PROJECT_ROOT>/processed/cellchat/smi_{P_T,P_TE,MP_T,MP_TE}_cellchat202511_updatecell2.0.qs
#             <PROJECT_ROOT>/processed/cellchat/smi_{P_T,P_TE,MP_T,MP_TE}_cellchat_copy202511_updatecell2.0.qs
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 01_cellchat_build_final_level4.R
# -----------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Dependencies
## ---------------------------------------------------------------------------
library(Seurat)
library(CellChat)
library(future)
library(qs)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
dir.create(file.path(PROJECT_ROOT, "processed", "cellchat"), recursive = TRUE, showWarnings = FALSE)

## Load the annotated object produced by 00_setup_and_annotation.R
smi <- qread(file.path(PROJECT_ROOT, "processed", "smi_annotated.qs"))


## ---------------------------------------------------------------------------
## Signalling database: keep only secreted signalling
## ---------------------------------------------------------------------------
CellChatDB <- CellChatDB.human  # use CellChatDB.mouse if running on mouse
# Subset CellChatDB for the cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB,
                           search = "Secreted Signaling",  # Cell-Cell Contact, ECM-Receptor, Non-protein Signaling, Secreted Signaling
                           key = "annotation")


## ===========================================================================
## 1.1 P_Tumour
## ===========================================================================
# Extract the Seurat object
smi_P_T <- subset(smi, subset = group == "P" & tissue == "T" & final_Level4_celltype != "Mix")
# Reset the factor levels of final_Level4_celltype
smi_P_T$final_Level4_celltype <- droplevels(smi_P_T$final_Level4_celltype)
levels(smi_P_T$final_Level4_celltype)
# CellChat needs normalised data
# Extract coordinates (relative to the FOV, each FOV is 510 um)
coords <- smi_P_T@meta.data[, c("x_FOV_px", "y_FOV_px")]
rownames(coords) <- rownames(smi_P_T@meta.data)
# A tolerance (tol) can be set to the radius or to a custom value
cellchat_spatial_factors <- data.frame(
  ratio = 0.12,   # coordinates are in px and must be scaled; the slide unit is
                  # mm while the position inside a FOV is in px, 1 px = 120 nm
  tol   = 15      # used as the tolerance
)
# Build the CellChat object
smi_P_T_cellchat <- createCellChat(
  object = smi_P_T,
  group.by = "final_Level4_celltype",   # cell-type grouping
  assay = "SCT",                        # the default assay is RNA; SCT is used here
  datatype = "spatial",                 # tell CellChat this is spatial data
  coordinates = coords,
  do.sparse = TRUE,                     # sparse matrix for speed
  spatial.factors = cellchat_spatial_factors)
## Set the ligand-receptor interaction database
smi_P_T_cellchat@DB <- CellChatDB.use
## Preprocess the expression data for communication analysis
# Subset the expression data of the signalling genes
smi_P_T_cellchat <- subsetData(smi_P_T_cellchat)
future::plan("multisession", workers = 40)  # multithreading
options(future.globals.maxSize = 100000 * 1024^5)
# Identify over-expressed genes
smi_P_T_cellchat <- identifyOverExpressedGenes(smi_P_T_cellchat)
# Identify over-expressed ligand-receptor pairs
smi_P_T_cellchat <- identifyOverExpressedInteractions(smi_P_T_cellchat)
# Project the gene expression data onto the PPI
# smi_P_T_cellchat <- smoothData(smi_P_T_cellchat, adj = PPI.human)  # smoothing may
#   reduce local spatial differences, so consider whether the raw spatial
#   features should be preserved

# No distance limit is used here; if there are no significant results, build one
# CellChat per sample and merge them for the computation.
smi_P_T_cellchat_copy <- smi_P_T_cellchat  # keep a copy
smi_P_T_cellchat <- computeCommunProb(
  smi_P_T_cellchat,
  type = "truncatedMean",  # mean gene expression per cell group; the default
                           # "triMean" yields fewer but stronger interactions,
                           # while "truncatedMean" needs a trim value and yields
                           # more interactions
  trim = 0.1, distance.use = TRUE, interaction.range = 250, scale.distance = 4,
  contact.dependent = TRUE, contact.range = 50)  # spatial distance constrains the communication probability
# Filter out cell-cell communication when a cell group has too few cells
smi_P_T_cellchat <- filterCommunication(smi_P_T_cellchat, min.cells = 10)
# Infer cell-cell communication at the signalling-pathway level
smi_P_T_cellchat <- computeCommunProbPathway(smi_P_T_cellchat)
# Compute the aggregated cell-cell communication network
smi_P_T_cellchat <- aggregateNet(smi_P_T_cellchat)
# Compute the network centrality scores
smi_P_T_cellchat <- netAnalysis_computeCentrality(
  smi_P_T_cellchat, slot.name = "netP")  # "netP" is the inferred pathway-level
                                         # communication network; the seed keeps results reproducible

# Save the data
qsave(smi_P_T_cellchat,      file.path(PROJECT_ROOT, "processed", "cellchat", "smi_P_T_cellchat202511_updatecell2.0.qs"))
qsave(smi_P_T_cellchat_copy, file.path(PROJECT_ROOT, "processed", "cellchat", "smi_P_T_cellchat_copy202511_updatecell2.0.qs"))
# smi_P_T_cellchat <- qread(file.path(PROJECT_ROOT, "processed", "cellchat", "smi_P_T_cellchat202511_updatecell2.0.qs"))


## ===========================================================================
## 1.2 P_TE (tumour margin)
## ===========================================================================
# Extract the Seurat object
smi_P_TE <- subset(smi, subset = group == "P" & tissue == "TE" & final_Level4_celltype != "Mix")
# Reset the factor levels of final_Level4_celltype
smi_P_TE$final_Level4_celltype <- droplevels(smi_P_TE$final_Level4_celltype)
levels(smi_P_TE$final_Level4_celltype)
# CellChat needs normalised data
# Extract coordinates (relative to the FOV, each FOV is 510 um)
coords <- smi_P_TE@meta.data[, c("x_FOV_px", "y_FOV_px")]
rownames(coords) <- rownames(smi_P_TE@meta.data)
# A tolerance (tol) can be set to the radius or to a custom value
cellchat_spatial_factors <- data.frame(
  ratio = 0.12,   # coordinates are in px and must be scaled; the slide unit is
                  # mm while the position inside a FOV is in px, 1 px = 120 nm
  tol   = 15      # used as the tolerance
)
# Build the CellChat object
smi_P_TE_cellchat <- createCellChat(
  object = smi_P_TE,
  group.by = "final_Level4_celltype",   # cell-type grouping
  assay = "SCT",                        # the default assay is RNA; SCT is used here
  datatype = "spatial",                 # tell CellChat this is spatial data
  coordinates = coords,
  do.sparse = TRUE,                     # sparse matrix for speed
  spatial.factors = cellchat_spatial_factors)
## Set the ligand-receptor interaction database
smi_P_TE_cellchat@DB <- CellChatDB.use
## Preprocess the expression data for communication analysis
# Subset the expression data of the signalling genes
smi_P_TE_cellchat <- subsetData(smi_P_TE_cellchat)
future::plan("multisession", workers = 40)  # multithreading
options(future.globals.maxSize = 100000 * 1024^5)
# Identify over-expressed genes
smi_P_TE_cellchat <- identifyOverExpressedGenes(smi_P_TE_cellchat)
# Identify over-expressed ligand-receptor pairs
smi_P_TE_cellchat <- identifyOverExpressedInteractions(smi_P_TE_cellchat)
# Project the gene expression data onto the PPI
# smi_P_TE_cellchat <- smoothData(smi_P_TE_cellchat, adj = PPI.human)  # smoothing may
#   reduce local spatial differences, so consider whether the raw spatial
#   features should be preserved

# No distance limit is used here; if there are no significant results, build one
# CellChat per sample and merge them for the computation.
smi_P_TE_cellchat_copy <- smi_P_TE_cellchat  # keep a copy
smi_P_TE_cellchat <- computeCommunProb(
  smi_P_TE_cellchat,
  type = "truncatedMean",  # mean gene expression per cell group; the default
                           # "triMean" yields fewer but stronger interactions,
                           # while "truncatedMean" needs a trim value and yields
                           # more interactions
  trim = 0.1, distance.use = TRUE, interaction.range = 200, scale.distance = 4,
  contact.dependent = TRUE, contact.range = 30)  # spatial distance constrains the communication probability
# Filter out cell-cell communication when a cell group has too few cells
smi_P_TE_cellchat <- filterCommunication(smi_P_TE_cellchat, min.cells = 10)
# Infer cell-cell communication at the signalling-pathway level
smi_P_TE_cellchat <- computeCommunProbPathway(smi_P_TE_cellchat)
# Compute the aggregated cell-cell communication network
smi_P_TE_cellchat <- aggregateNet(smi_P_TE_cellchat)
# Compute the network centrality scores
smi_P_TE_cellchat <- netAnalysis_computeCentrality(
  smi_P_TE_cellchat, slot.name = "netP")

# Save the data
qsave(smi_P_TE_cellchat,      file.path(PROJECT_ROOT, "processed", "cellchat", "smi_P_TE_cellchat202511_updatecell2.0.qs"))
qsave(smi_P_TE_cellchat_copy, file.path(PROJECT_ROOT, "processed", "cellchat", "smi_P_TE_cellchat_copy202511_updatecell2.0.qs"))
# smi_P_TE_cellchat <- qread(file.path(PROJECT_ROOT, "processed", "cellchat", "smi_P_TE_cellchat202511_updatecell2.0.qs"))


## ===========================================================================
## 1.3 MP_Tumour
## ===========================================================================
# Extract the Seurat object
smi_MP_T <- subset(smi, subset = group == "MP" & tissue == "T" & final_Level4_celltype != "Mix")
# Reset the factor levels of final_Level4_celltype
smi_MP_T$final_Level4_celltype <- droplevels(smi_MP_T$final_Level4_celltype)
levels(smi_MP_T$final_Level4_celltype)
# CellChat needs normalised data
# Extract coordinates (relative to the FOV, each FOV is 510 um)
coords <- smi_MP_T@meta.data[, c("x_FOV_px", "y_FOV_px")]
rownames(coords) <- rownames(smi_MP_T@meta.data)
# A tolerance (tol) can be set to the radius or to a custom value
cellchat_spatial_factors <- data.frame(
  ratio = 0.12,   # coordinates are in px and must be scaled; the slide unit is
                  # mm while the position inside a FOV is in px, 1 px = 120 nm
  tol   = 15      # used as the tolerance
)
# Build the CellChat object
smi_MP_T_cellchat <- createCellChat(
  object = smi_MP_T,
  group.by = "final_Level4_celltype",   # cell-type grouping
  assay = "SCT",                        # the default assay is RNA; SCT is used here
  datatype = "spatial",                 # tell CellChat this is spatial data
  coordinates = coords,
  do.sparse = TRUE,                     # sparse matrix for speed
  spatial.factors = cellchat_spatial_factors)
## Set the ligand-receptor interaction database
smi_MP_T_cellchat@DB <- CellChatDB.use
## Preprocess the expression data for communication analysis
# Subset the expression data of the signalling genes
smi_MP_T_cellchat <- subsetData(smi_MP_T_cellchat)
future::plan("multisession", workers = 40)  # multithreading
options(future.globals.maxSize = 100000 * 1024^5)
# Identify over-expressed genes
smi_MP_T_cellchat <- identifyOverExpressedGenes(smi_MP_T_cellchat)
# Identify over-expressed ligand-receptor pairs
smi_MP_T_cellchat <- identifyOverExpressedInteractions(smi_MP_T_cellchat)
# Project the gene expression data onto the PPI
# smi_MP_T_cellchat <- smoothData(smi_MP_T_cellchat, adj = PPI.human)  # smoothing may
#   reduce local spatial differences, so consider whether the raw spatial
#   features should be preserved

# No distance limit is used here; if there are no significant results, build one
# CellChat per sample and merge them for the computation.
smi_MP_T_cellchat_copy <- smi_MP_T_cellchat  # keep a copy
smi_MP_T_cellchat <- computeCommunProb(
  smi_MP_T_cellchat,
  type = "truncatedMean",  # mean gene expression per cell group; the default
                           # "triMean" yields fewer but stronger interactions,
                           # while "truncatedMean" needs a trim value and yields
                           # more interactions
  trim = 0.1, distance.use = TRUE, interaction.range = 250, scale.distance = 4,
  contact.dependent = TRUE, contact.range = 50)  # spatial distance constrains the communication probability
# Filter out cell-cell communication when a cell group has too few cells
smi_MP_T_cellchat <- filterCommunication(smi_MP_T_cellchat, min.cells = 10)
# Infer cell-cell communication at the signalling-pathway level
smi_MP_T_cellchat <- computeCommunProbPathway(smi_MP_T_cellchat)
# Compute the aggregated cell-cell communication network
smi_MP_T_cellchat <- aggregateNet(smi_MP_T_cellchat)
# Compute the network centrality scores
smi_MP_T_cellchat <- netAnalysis_computeCentrality(
  smi_MP_T_cellchat, slot.name = "netP")

# Save the data
qsave(smi_MP_T_cellchat,      file.path(PROJECT_ROOT, "processed", "cellchat", "smi_MP_T_cellchat202511_updatecell2.0.qs"))
qsave(smi_MP_T_cellchat_copy, file.path(PROJECT_ROOT, "processed", "cellchat", "smi_MP_T_cellchat_copy202511_updatecell2.0.qs"))
# smi_MP_T_cellchat <- qread(file.path(PROJECT_ROOT, "processed", "cellchat", "smi_MP_T_cellchat202511_updatecell2.0.qs"))


## ===========================================================================
## 1.4 MP_TE (tumour margin)
## ===========================================================================
# Extract the Seurat object
smi_MP_TE <- subset(smi, subset = group == "MP" & tissue == "TE" & final_Level4_celltype != "Mix")
# Reset the factor levels of final_Level4_celltype
smi_MP_TE$final_Level4_celltype <- droplevels(smi_MP_TE$final_Level4_celltype)
levels(smi_MP_TE$final_Level4_celltype)
# CellChat needs normalised data
# Extract coordinates (relative to the FOV, each FOV is 510 um)
coords <- smi_MP_TE@meta.data[, c("x_FOV_px", "y_FOV_px")]
rownames(coords) <- rownames(smi_MP_TE@meta.data)
# A tolerance (tol) can be set to the radius or to a custom value
cellchat_spatial_factors <- data.frame(
  ratio = 0.12,   # coordinates are in px and must be scaled; the slide unit is
                  # mm while the position inside a FOV is in px, 1 px = 120 nm
  tol   = 15      # used as the tolerance
)
# Build the CellChat object
smi_MP_TE_cellchat <- createCellChat(
  object = smi_MP_TE,
  group.by = "final_Level4_celltype",   # cell-type grouping
  assay = "SCT",                        # the default assay is RNA; SCT is used here
  datatype = "spatial",                 # tell CellChat this is spatial data
  coordinates = coords,
  do.sparse = TRUE,                     # sparse matrix for speed
  spatial.factors = cellchat_spatial_factors)
## Set the ligand-receptor interaction database
smi_MP_TE_cellchat@DB <- CellChatDB.use
## Preprocess the expression data for communication analysis
# Subset the expression data of the signalling genes
smi_MP_TE_cellchat <- subsetData(smi_MP_TE_cellchat)
future::plan("multisession", workers = 40)  # multithreading
options(future.globals.maxSize = 100000 * 1024^5)
# Identify over-expressed genes
smi_MP_TE_cellchat <- identifyOverExpressedGenes(smi_MP_TE_cellchat)
# Identify over-expressed ligand-receptor pairs
smi_MP_TE_cellchat <- identifyOverExpressedInteractions(smi_MP_TE_cellchat)
# Project the gene expression data onto the PPI
# smi_MP_TE_cellchat <- smoothData(smi_MP_TE_cellchat, adj = PPI.human)  # smoothing may
#   reduce local spatial differences, so consider whether the raw spatial
#   features should be preserved

# No distance limit is used here; if there are no significant results, build one
# CellChat per sample and merge them for the computation.
smi_MP_TE_cellchat_copy <- smi_MP_TE_cellchat  # keep a copy
smi_MP_TE_cellchat <- computeCommunProb(
  smi_MP_TE_cellchat,
  type = "truncatedMean",  # mean gene expression per cell group; the default
                           # "triMean" yields fewer but stronger interactions,
                           # while "truncatedMean" needs a trim value and yields
                           # more interactions
  trim = 0.1, distance.use = TRUE, interaction.range = 200, scale.distance = 4,
  contact.dependent = TRUE, contact.range = 30)  # spatial distance constrains the communication probability
# Filter out cell-cell communication when a cell group has too few cells
smi_MP_TE_cellchat <- filterCommunication(smi_MP_TE_cellchat, min.cells = 10)
# Infer cell-cell communication at the signalling-pathway level
smi_MP_TE_cellchat <- computeCommunProbPathway(smi_MP_TE_cellchat)
# Compute the aggregated cell-cell communication network
smi_MP_TE_cellchat <- aggregateNet(smi_MP_TE_cellchat)
# Compute the network centrality scores
smi_MP_TE_cellchat <- netAnalysis_computeCentrality(
  smi_MP_TE_cellchat, slot.name = "netP")

# Save the data
qsave(smi_MP_TE_cellchat,      file.path(PROJECT_ROOT, "processed", "cellchat", "smi_MP_TE_cellchat202511_updatecell2.0.qs"))
qsave(smi_MP_TE_cellchat_copy, file.path(PROJECT_ROOT, "processed", "cellchat", "smi_MP_TE_cellchat_copy202511_updatecell2.0.qs"))
# smi_MP_TE_cellchat <- qread(file.path(PROJECT_ROOT, "processed", "cellchat", "smi_MP_TE_cellchat202511_updatecell2.0.qs"))
