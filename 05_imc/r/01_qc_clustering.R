# -----------------------------------------------------------------------------
# QC, clustering and dimensionality reduction of IMC single-cell data
#
# Purpose   : Build a SpatialExperiment from the per-ROI segmentation CSVs,
#             apply the arcsinh / quantile / min-max normalisation used
#             throughout the study, run SVD-based PCA + Harmony batch
#             correction, UMAP and PhenoGraph (FastPG) clustering. This is the
#             upstream skeleton for every downstream IMC analysis.
#
# Inputs    : <PROJECT_ROOT>/imc/csv/*.csv        per-cell segmentation tables
#             <PROJECT_ROOT>/imc/info.csv         sample / ROI / group metadata
#             <PROJECT_ROOT>/imc/marker.csv       antibody panel (channel, marker, lineage)
#
# Outputs   : <PROJECT_ROOT>/imc/out/0_QC/CellCountsPerROI.csv
#             <PROJECT_ROOT>/imc/out/1_Cluster/1_spe.rds
#             <PROJECT_ROOT>/imc/out/1_Cluster/annotation.csv
#             <PROJECT_ROOT>/imc/out/1_Cluster/1.3_Heatmap/*.csv
#             <PROJECT_ROOT>/imc/out/1_Cluster/1.5_echo_roi/freq.csv
#             <PROJECT_ROOT>/imc/out/1_Cluster/1.6_roi.anno/*.csv
#             <PROJECT_ROOT>/imc/out/1_Cluster/1.7_MarkerNames/*.csv
#             <PROJECT_ROOT>/imc/out/1_Cluster/1.8_RoiMarker/*.csv
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 01_qc_clustering.R
# -----------------------------------------------------------------------------

library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)
library(ggpubr)
library(plyr)
library(ggbeeswarm)
library(scales)
library(cowplot)
library(viridis)
library(scater)
library(BiocParallel)
library(harmony)
library(SingleCellExperiment)
library(SpatialExperiment)
library(tidyverse)
library(imcRtools)
library(diffcyt)
library(dittoSeq)
library(bluster)
library(argparser)
library(FastPG)
library(glue)
library(fs)
set.seed(100)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

csv.dir <- file.path(PROJECT_ROOT, "imc", "csv")
info.file <- file.path(PROJECT_ROOT, "imc", "info.csv")
marker.file <- file.path(PROJECT_ROOT, "imc", "marker.csv")
output.dir <- file.path(PROJECT_ROOT, "imc", "out")
if(! dir.exists(output.dir)) {dir.create(output.dir)}


group.stat.list <- list(c("T","N"),c("T","MT"),c("N","MT"))

cofactor = 1
q.nor.n = 101




## 0. QC ##########################################################
output.dir0 <-  paste(output.dir, '0_QC', sep = '/')
if(! dir.exists(output.dir0)) {dir.create(output.dir0)}

## read marker panel
panel <- read.csv(marker.file)
colnames(panel) <- c('channel', 'marker','lineager')


# read csv
files <- list.files(csv.dir ,pattern='.csv$', full=TRUE)
csv.list <- lapply(files, function(x) {
  temp <- read.csv(x)
  col.temp <- colnames(temp)[grep('X[0-9]+', colnames(temp))]
  col.channel <- str_split( col.temp , pattern = '_', simplify = T)[,1] %>% str_sub(2, nchar(.))
  colnames(temp)[grep('X[0-9]+', colnames(temp))]  <- panel$marker[match(col.channel, panel$channel)]
  temp <- temp[, !is.na(colnames(temp))]
  temp
})
names(csv.list) <- files

## check file colnames
for(x in files) {
  temp <- csv.list[[x]] %>% colnames()
  temp <- any(str_detect(temp,'[.]1$'))
  if(temp == TRUE) {stop( c(x,'colnames are uncorrect'))}
  if(is.na(temp)) {stop( c('Panel marker file is uncorrect'))}
}
cur_features <- do.call('rbind.fill', csv.list)


## read sample info ####
info <- read.csv(info.file)
colnames(info)[1] <- 'sample_id'

info<- info[,-grep("order", names(info)),drop=F]

if (!'roi_location' %in% (info %>% colnames())){
  filess <- list.files(csv.dir ,pattern='.csv$')

  TName <- str_extract_all(filess,pattern = 'T[0-9]+', simplify = T)
  ROIName <- str_extract_all(filess,pattern = 'ROI_[A-Z0-9_A-Z0-9]+', simplify = T)
  ROIID <- str_extract_all(filess,pattern = 'ROI_[A-Z0-9_A-Z0-9]+', simplify = T)

  matchdata<-data.frame(roi_id=paste0(TName,'_',str_replace(ROIID[,1],'_','')),
                        roi_location=paste0(TName,'_',str_replace(ROIName[,1],'_','')))
  info<-merge(info,matchdata)
}

TName <- sapply(csv.list, nrow) %>%
  data.frame() %>%  rownames() %>%
  str_extract_all(pattern = 'T[0-9]+', simplify = T)
ROIName <- sapply(csv.list, nrow) %>%
  data.frame() %>%  rownames() %>%
  str_extract_all(pattern = 'ROI_[A-Z0-9_A-Z0-9]+', simplify = T)
roi_id.csv <- paste(TName[,1], str_replace(ROIName[,1],'_',''), sep = '_')

temp_Tname <- str_extract_all(cur_features$roi, pattern = 'T[0-9]+', simplify = T)
temp_ROIname <- str_extract_all(cur_features$roi, pattern = 'ROI_[A-Z0-9_A-Z0-9]+', simplify = T)
temp <- paste(temp_Tname[,1], str_replace(temp_ROIname[,1],'_',''), sep = '_')
cur_features$roi <- temp


roi_info <- data.frame(number = sapply(csv.list, nrow), roi_id = roi_id.csv) %>%
  group_by(roi_id) %>%
  left_join(info) %>% mutate(CountsPermm2 = number / roi_area)

colnames(roi_info)[1]<-'counts'
roi_info %>% write.csv(paste(output.dir0, 'CellCountsPerROI.csv', sep = '/'))


row.idx <- is.na(roi_info) %>% apply(1, any) %>% which
if(length(row.idx) > 0) {stop(paste0( paste(roi_info$roi_id[row.idx], collapse = ' '),' info is missing.'))}


group.num <- colnames(roi_info)[grep('group',colnames(roi_info))]

group.num <- colnames(info)[grep('group',colnames(info))]

marker <- panel %>% filter(! lineager == 0) %>% pull(marker)
counts <- cur_features %>% select(marker)

meta <- cur_features %>% select(roi, AreaID,  CellID) %>%
  mutate(roi_id = roi, .keep = 'unused') %>%
  left_join(info) ### add sample info
coords <- cur_features %>% select(contains('position'))
colnames(coords) <- c("Pos_X", "Pos_Y")



spe <- SpatialExperiment(assays = list(counts = t(counts)),
                         colData = meta,
                         sample_id = as.character(meta$sample_id),
                         image_id = as.character(meta$roi_id),
                         spatialCoords = as.matrix(coords))
colnames(spe) <- paste0(spe$roi_id, ".", spe$CellID)


rowData(spe)$use_channel <- rownames(spe) %in% (panel %>% filter(lineager == '1') %>% pull(marker))

counts.tsf <- asinh(counts(spe)/cofactor)
counts.nor <- counts.tsf
for (i in 1:nrow(counts.nor)) {
  print(i)
  counts.nor[i,] <- q.nor(x = counts.nor[i,], n = q.nor.n)
  counts.nor[i,] <- minMax(counts.nor[i,])
}
assay(spe, "exprs") <- counts.nor

set.seed(220225)

spe <- runUMAP(spe, subset_row = rowData(spe)$use_channel, exprs_values = "exprs")
mat <- t(assay(spe, "exprs")) [,rowData(spe)$use_channel]
svd_res <- svd(mat)
variance_explained <- (svd_res$d^2)/sum(svd_res$d^2)
npc <- which(cumsum(variance_explained) >= 0.95)[1]


harmony_emb <- HarmonyMatrix(mat, as.factor(spe$batch_id), do_pca = T, npcs = npc) ### batch_id
reducedDim(spe, "harmony") <- harmony_emb
spe <- runUMAP(spe, dimred = "harmony", name = "UMAP_harmony")


# plot
spe$batch_id <- factor(spe$batch_id, levels = unique(spe$batch_id))

temp <- levels(spe$batch_id) %>% length
ncol.n <- ceiling(temp/10)


 ## 1. reduction and cluster #############################################################
print('start section 1: reduction and cluster')
output.dir1 <-  paste(output.dir, '1_Cluster', sep = '/')
if(! dir.exists(output.dir1)) {dir.create(output.dir1)}

### 1.1 cluster ####
expr.mat <- reducedDim(spe, "harmony")
cluster_PhenoGraph <- FastPG::fastCluster(as.matrix(expr.mat),k = 30, num_threads = 100)
spe$pg_cluster <- as.factor(cluster_PhenoGraph$communities)

temp <- data.frame(spe$pg_cluster, reducedDim(spe, "UMAP_harmony"))
colnames(temp) <- c('pg_cluster','UMAP_1','UMAP_2')
lc.cent = temp %>% group_by(pg_cluster) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarise_all(median)

### 1.2 feature plot ####
output.dir1sub <-  paste(output.dir1, '1.2_FeaturePlot', sep = '/')
if(! dir.exists(output.dir1sub)) {dir.create(output.dir1sub)}

### 1.3 marker median heatmap ####
hmdata <- cbind(t(assay(spe, "exprs"))[,rowData(spe)$use_channel], pg_cluster = spe$pg_cluster) %>%
  as_tibble() %>% group_by(pg_cluster) %>%
  summarise_all(mean) %>% column_to_rownames(var = 'pg_cluster')

rownames(hmdata) <- spe$pg_cluster %>% levels()


dir.create(paste0(output.dir1, '/1.3_Heatmap'))
hmdata %>% write.csv(paste(output.dir1, '/1.3_Heatmap/heatmap.median.csv', sep = '/'))

hmdata <- cbind(t(assay(spe, "exprs")), pg_cluster = spe$pg_cluster) %>%
  as_tibble() %>% group_by(pg_cluster) %>%
  summarise_all(mean) %>% column_to_rownames(var = 'pg_cluster')

rownames(hmdata) <- spe$pg_cluster %>% levels()

hmdata %>% write.csv(paste(output.dir1, '/1.3_Heatmap/heatmap.median.all.csv', sep = '/'))


celllabel<-data.frame(pg_cluster=levels(spe$pg_cluster),
                      new_pg_cluster=rep('',length(levels(spe$pg_cluster))),
                      celltype=rep('',length(levels(spe$pg_cluster))))
write.csv(celllabel,file = paste0(output.dir1,'/annotation.csv'))
spe@metadata<-grouporder
saveRDS(spe, paste(output.dir1, "1_spe.rds", sep = '/'))


outdir11 <- paste0(output.dir1,'/1.5_echo_roi')
if(!dir.exists(outdir11)) {dir.create(outdir11)}

outdir22 <- paste0(output.dir1,'/1.6_roi.anno')
if(!dir.exists(outdir22)) {dir.create(outdir22)}

## freq table
freqROI <- colData(spe) %>% as_tibble() %>% dplyr::select(roi_location, pg_cluster)
freqROI$pg_cluster <- paste0('cluster',freqROI$pg_cluster)
freqROI%>% group_by(roi_location) %>% table() %>%
  write.csv(paste0(outdir11,'/freq.csv'))

## roi anno
temp <- cbind(spe@colData,spe@int_colData$spatialCoords) %>% as_tibble() %>% select(roi_location, CellID, pg_cluster,Pos_X,Pos_Y)
colnames(temp)[3] <- 'CellType'
temp$CellType <- paste0('cluster',temp$CellType)
temp.list <- split(temp, temp$roi_location)


lapply(temp.list, function(x) {
  x %>% dplyr::select(CellID, CellType,Pos_X,Pos_Y) %>%
    arrange(CellID) %>%
    write.csv(paste(outdir22, paste0(unique(x$roi_location),'.csv'), sep = '/'), row.names = F)
})

marker.df <- t(assay(spe, "exprs")) %>% as.tibble() %>% mutate(roi_location = spe$roi_location,CellID=spe$CellID,Pos_X=spe@int_colData$spatialCoords[,1],Pos_Y=spe@int_colData$spatialCoords[,2])

MarkerName <- colnames(marker.df)
MarkerName <- MarkerName[!MarkerName %in% c("roi_location","CellID")]

dir_create(paste(output.dir1,'1.7_MarkerNames',sep = '/'))
MarkerName %>% write.csv(paste(output.dir1,'1.7_MarkerNames','MarkerNames.csv',sep = '/'),row.names = F)

roiList <- marker.df$roi_location %>% unique()
dir_create(paste(output.dir1,'1.8_RoiMarker',sep = '/'))
for (roi in roiList){
  marker.df.sub <- marker.df[marker.df$roi_location %in% roi,]
  marker.df.sub %>% write.csv(paste(output.dir1,'1.8_RoiMarker',glue('{roi}.csv'),sep = '/'),row.names = F)

}
print('QC_Cluster Done' )
