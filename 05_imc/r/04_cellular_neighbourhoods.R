# -----------------------------------------------------------------------------
# Cellular neighbourhoods (CN) by k-means on aggregated neighbour phenotypes
#
# Purpose   : Build a k-NN spatial graph, aggregate the cell-type composition
#             of the neighbours of each cell and cluster those neighbourhood
#             profiles with k-means to define cellular neighbourhoods (CN). CN
#             composition, per-ROI CN frequencies and the group statistics are
#             written for the 10 / 15 / 20 cluster solutions.
#
# Inputs    : <PROJECT_ROOT>/imc/out/1_Cluster/1_spe.rds
#
# Outputs   : <PROJECT_ROOT>/imc/out/4_CN/CN_<n>.rds
#             <PROJECT_ROOT>/imc/out/4_CN/CN_<n>/4.3_CN_Heatmap_<n>/*.csv
#             <PROJECT_ROOT>/imc/out/4_CN/CN_<n>/4.2_CN_ratio_<n>/CN_ratio.csv
#             <PROJECT_ROOT>/imc/out/4_CN/CN_<n>/4.2_CN_ratio_<n>/CN_count.csv
#             <PROJECT_ROOT>/imc/out/4_CN/CN_<n>/4.2_CN_ratio_<n>/<group>/*.csv
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 04_cellular_neighbourhoods.R
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
library(FastPG)
library(rstatix)
library(paletteer)
library(fs)
library(Seurat)
library(glue)
library(ggvoronoi)

set.seed(100)
palette <- c(paletteer_d("ggsci::category20_d3"),paletteer_d("ggsci::category20b_d3"),paletteer_d("ggsci::category20c_d3"))

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

rds <- file.path(PROJECT_ROOT, "imc", "out", "1_Cluster", "1_spe.rds")
spe <- readRDS(rds)
grouporder<-spe@metadata
graphType <- "knn"

# output dir
output.dir <- file.path(PROJECT_ROOT, "imc", "out")
if(! dir.exists(output.dir)) {dir.create(output.dir)}

group.stat.list <- list(c("T","N"),c("T","MT"),c("N","MT"))
## stat
knn_k = 20



## 3. CN #########################################################################
output.dir3 <-  paste(output.dir, '4_CN', sep = '/')
if(! dir.exists(output.dir3)) {dir.create(output.dir3)}
set.seed(220705)


if (graphType=='knn'){
  spe <- buildSpatialGraph(spe, img_id = "roi_id", type = "knn", k = knn_k)
  spe <- aggregateNeighbors(spe, colPairName = "knn_interaction_graph",
                          aggregate_by = "metadata", count_by = "celltype")
}else if (graphType=='expansion'){
  spe <- buildSpatialGraph(spe, img_id = "roi_id", type = "expansion", threshold = knn_k)
  spe <- aggregateNeighbors(spe, colPairName = "expansion_interaction_graph",
                            aggregate_by = "metadata", count_by = "celltype")
}


roi.num <- spe$roi_id %>% unique() %>% length

## CN number #####
for(CN.num in c(10,15,20)) {
  print(CN.num)
  # CN.num = as.numeric(CN.num)

  dir.create(paste0(output.dir3,'/CN_',CN.num))

  dir_create(paste0(output.dir3,'/CN_',CN.num,'/4.1_CN_graph','_',CN.num))
  dir_create(paste0(output.dir3,'/CN_',CN.num,'/4.2_CN_ratio','_',CN.num))
  dir.create(paste0(output.dir3,'/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num))

  output.dir4 <-  paste(output.dir3, paste0('/CN_',CN.num,'/4.2_CN_ratio','_',CN.num), sep = '/') ## modified


  cn_1 <- kmeans(spe$aggregatedNeighbors, centers = CN.num)
  #cn_1 <- kmeans(spe$aggregatedNeighbors, centers = 1)
  spe$cn_celltypes <- as.factor(cn_1$cluster)

  saveRDS(spe, file = paste0(output.dir3,'/CN_',CN.num,'.rds'))  # *

  if (is.null(spe$roi_location)){
    spe$roi_location <- spe$roi_id
  }

  for_plot <- colData(spe) %>% as_tibble() %>%
    group_by(cn_celltypes, celltype) %>%
    dplyr::summarize(count = n()) %>%
    mutate(freq = count / sum(count)) %>%
    pivot_wider(id_cols = cn_celltypes, names_from = celltype,
                values_from = freq, values_fill = 0) %>%
    ungroup() %>%
    column_to_rownames('cn_celltypes')
  rownames(for_plot) <- paste0('CN', rownames(for_plot))
  for_plot %>% write.csv(paste(output.dir3, paste0('/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num,'/','CN_heatmap_',CN.num,'_rowscale.csv'), sep = '/'))


  #paste0(output.dir3,'/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num)

  for_plot <- colData(spe) %>% as_tibble() %>%
    group_by( celltype,cn_celltypes) %>%
    dplyr::summarize(count = n()) %>%
    mutate(freq = count / sum(count)) %>%
    pivot_wider(id_cols = cn_celltypes, names_from = celltype,
                values_from = freq, values_fill = 0) %>%
    ungroup() %>%
    column_to_rownames('cn_celltypes')
  rownames(for_plot) <- paste0('CN', rownames(for_plot))
  for_plot %>% write.csv(paste(output.dir3, paste0('/CN_',CN.num,'/4.3_CN_Heatmap','_',CN.num,'/','CN_heatmap_',CN.num,'_colscale.csv'), sep = '/'))



  #output.dir4 <-  paste(output.dir3, paste0('CN.ratio_',1), sep = '/')
  if(! dir.exists(output.dir4)) {dir.create(output.dir4)}

  meta.info <- colData(spe) %>% as_tibble() %>% select(roi_id, contains('group')) %>% unique()
  ttest.data <- colData(spe) %>% as_tibble() %>% group_by(roi_id) %>% dplyr::select(cn_celltypes)%>%
    table()%>% normalize()%>% as.data.frame() %>% left_join(meta.info)
  ttest.data %>% write.csv(file = paste(output.dir4, 'CN_ratio.csv', sep = '/'))

  cellcount<-colData(spe) %>% as_tibble() %>% group_by(roi_id) %>% dplyr::select(cn_celltypes)%>%
    table()%>% as.data.frame() %>% left_join(meta.info)
   names(cellcount)[3]<-'Count'
    write.csv(cellcount,file = paste(output.dir4, 'CN_count.csv', sep = '/'))




  ### 4.1 t test #####
  group.num <- colnames(ttest.data)[grep('group', colnames(ttest.data))]
  for (g in group.num) {

    stat.dir <- paste(output.dir4, g, sep = '/')
    if(! dir.exists(stat.dir)) {dir.create(stat.dir)}

    ######group in CN
    histgram.data <- colData(spe) %>% as.data.frame() %>% select(cn_celltypes,g)
    colnames(histgram.data)[2] <- 'group'
    histgram.data<-subset(histgram.data,!group %in% '')
    histgram.data<-table(histgram.data$group,histgram.data$cn_celltypes) %>% prop.table(margin = 2) %>% as.data.frame()
    names(histgram.data)[1:2]<-c('Group','CN')
    histgram.data$CN<-paste0('CN',histgram.data$CN)
    histgram.data$CN<-factor(histgram.data$CN,levels = unique(histgram.data$CN))


    write.csv(histgram.data,file = paste(output.dir4, g,'Ratio_group_in_CN.csv', sep = '/'))

    histgram.data <- colData(spe) %>% as.data.frame() %>% select(celltype,cn_celltypes)
    colnames(histgram.data)[2] <- 'CN'
    histgram.data<-table(histgram.data$celltype,histgram.data$CN) %>% prop.table(margin = 2) %>% as.data.frame()
    names(histgram.data)[1:2]<-c('Celltype','CN')
    histgram.data$CN<-paste0('CN',histgram.data$CN)
    histgram.data$CN<-factor(histgram.data$CN,levels = unique(histgram.data$CN))


    write.csv(histgram.data,file = paste(output.dir4, g,'Ratio_celltype_in_CN.csv', sep = '/'))




    ##### celltype in group in CN
    histgram.data <- colData(spe) %>% as.data.frame() %>% select(cn_celltypes,g,celltype)
    colnames(histgram.data)[1:2] <- c('CN','group')
    histgram.data<-subset(histgram.data,!group %in% '')
    frequency_table <- histgram.data %>%
      group_by(CN, group, celltype) %>%
      summarise(Frequency = n(), .groups = 'drop')

    # compute relative frequencies
    frequency_table <- frequency_table %>%
      group_by(CN, group) %>%
      mutate(RelativeFrequency = Frequency / sum(Frequency)) %>%
      ungroup()

    write.csv(frequency_table,file = paste(output.dir4, g,'Ratio_ct _in_group_in_CN.csv', sep = '/'))

  }
}
print(':: CN has already done ::')
