# -----------------------------------------------------------------------------
# Cell-type proportions per ROI
#
# Purpose   : Compute the cell-type composition (relative frequency) and the
#             raw cell counts of every ROI from the annotated SpatialExperiment
#             and write the per-ROI tables used for the region-wise frequency
#             plots.
#
# Inputs    : <PROJECT_ROOT>/imc/out/1_Cluster/1_spe.rds
#
# Outputs   : <PROJECT_ROOT>/imc/out/3_Cluster.ratio/<category>/cluster_ratio.csv
#             <PROJECT_ROOT>/imc/out/3_Cluster.ratio/<category>/cluster_count.csv
#             <PROJECT_ROOT>/imc/out/3_Cluster.ratio/<category>/3_spe_<category>.rds
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 03_cell_type_proportions.R
# -----------------------------------------------------------------------------

library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)
library(ggpubr)
library(plyr)
library(ggbeeswarm)
library(argparser)
library(scales)
library(SingleCellExperiment)
library(SpatialExperiment)
library(imcRtools)
library(diffcyt)
library(cowplot)
library(dittoSeq)
library(viridis)
library(bluster)
library(BiocParallel)
library(tidyverse)
library(scater)
library(harmony)
library(stringr)
# library(FastPG)
library(glue)
library(fs)
set.seed(100)

normalize <- function(x) {x/rowSums(x)}
z.score <- function(x) {(x-mean(x))/sd(x)}

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

# output dir
output.dir <- file.path(PROJECT_ROOT, "imc", "out")
if(! dir.exists(output.dir)) {dir.create(output.dir)}

output.dir1 <- paste0(output.dir, '/3_Cluster.ratio')
if(! dir.exists(output.dir1)) {dir.create(output.dir1)}

rds <- file.path(PROJECT_ROOT, "imc", "out", "1_Cluster", "1_spe.rds")
spe <- readRDS(rds)
grouporder<-spe@metadata
ori.all <- spe  # save original spe

## stat

group.stat.list <- list(c("T","N"),c("T","MT"),c("N","MT"))

if (is.null(spe$roi_location)){
  spe$roi_location <- spe$roi_id
}




partition_by <- argv$partition_by
if(partition_by %in% colnames(spe@colData)) {
  prefix <- spe@colData[[partition_by]] %>% unique()
} else {
  prefix <- 'all'
}


for (category in prefix) {
  print(category)
  if(category != 'all') {spe <- ori.all[,ori.all@colData$Imm_or_not %in% category]}

  output.dir2 <- paste0(output.dir1,'/',category)
  if(! dir.exists(output.dir2)) {dir.create(output.dir2)}




## 3. cell ratio stat #######################################################################

meta.info <- colData(spe) %>% as_tibble() %>% select(roi_location,roi_id, contains('group')) %>% unique()

ttest.data <- colData(spe) %>% as_tibble() %>% group_by(roi_id) %>% dplyr::select(celltype)%>%
  table() %>% normalize()%>% as.data.frame() %>% left_join(meta.info)

ttest.data %>% write.csv(file = paste(output.dir2, 'cluster_ratio.csv', sep = '/'))

colData(spe) %>% as_tibble() %>% group_by(roi_id) %>% dplyr::select(celltype)%>%
  table()%>% as.data.frame() %>% left_join(meta.info) %>%
  write.csv(file = paste(output.dir2, 'cluster_count.csv', sep = '/'))


saveRDS(spe, paste0(output.dir2, "/3_spe_", category, ".rds"))}
print('all done.')
