# -----------------------------------------------------------------------------
# Cell-type annotation of IMC clusters and lineage marker heatmaps
#
# Purpose   : Map the PhenoGraph clusters (from 01_qc_clustering.R) to named
#             cell lineages using a manually curated annotation.csv, derive the
#             annotated UMAP, the lineage-level marker heatmaps and the per-ROI
#             cell-type tables used for all downstream spatial analyses.
#
# Inputs    : <PROJECT_ROOT>/imc/out/1_Cluster/1_spe.rds   SpatialExperiment from step 01
#             <PROJECT_ROOT>/imc/annotation.csv            cluster -> cell type mapping
#             <PROJECT_ROOT>/imc/marker-f.csv              functional marker panel
#
# Outputs   : <PROJECT_ROOT>/imc/out/2_CellType/2_spe.rds
#             <PROJECT_ROOT>/imc/out/2_CellType/2.2_Heatmap/*.csv
#             <PROJECT_ROOT>/imc/out/2_CellType/2.8_cluster.marker/<group>/expresion.csv
#             <PROJECT_ROOT>/imc/out/2_CellType/2.4_echo_roi/freq.csv
#             <PROJECT_ROOT>/imc/out/2_CellType/2.5_roi.anno/*.csv
#             <PROJECT_ROOT>/imc/out/2_CellType/2.6_MarkerNames/*.csv
#             <PROJECT_ROOT>/imc/out/2_CellType/2.7_RoiMarker/*.csv
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 02_cell_type_annotation.R
# -----------------------------------------------------------------------------

library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)
library(ggpubr)
library(plyr)
library(tibble)
library(rstatix)
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
#library(FastPG)
library(glue)
library(fs)
set.seed(100)

normalize <- function(x) {x/rowSums(x)}
z.score <- function(x) {(x-mean(x))/sd(x)}


PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

anno.file <- file.path(PROJECT_ROOT, "imc", "annotation.csv")
FunctionMarkerCsvDir <- marker.file <- file.path(PROJECT_ROOT, "imc", "marker-f.csv")
rds <- file.path(PROJECT_ROOT, "imc", "out", "1_Cluster", "1_spe.rds")


spe <- readRDS(rds)

grouporder<-spe@metadata
output.dir <- file.path(PROJECT_ROOT, "imc", "out")
if(! dir.exists(output.dir)) {dir.create(output.dir)}

group.stat.list <- list(c("T","N"),c("T","MT"),c("N","MT"))


###### 3.1 Annotation #############################################################
output.dir1 <-  paste(output.dir, '2_CellType', sep = '/')
if(! dir.exists(output.dir1)) {dir.create(output.dir1)}

celltype <- read.csv(anno.file)
print(celltype)

spe$celltype<-mapvalues(spe$pg_cluster,from = celltype$new_pg_cluster,to = celltype$celltype)
spe <- spe[, ! spe$celltype == 'rm']
spe$celltype<-factor(spe$celltype,levels = unique(celltype$celltype)[which(unique(celltype$celltype)!='rm')])

spe$new_pg_cluster<-mapvalues(spe$pg_cluster,from = celltype$pg_cluster,to = celltype$new_pg_cluster)
celltype_number<-unique(celltype$celltype)[which(unique(celltype$celltype)!='rm')] %>% length
print(celltype_number)

celltype_umap<-paste0(c(1:celltype_number),':',unique(celltype$celltype)[which(unique(celltype$celltype)!='rm')])
spe$celltype_umap<-mapvalues(spe$celltype,from = unique(celltype$celltype)[which(unique(celltype$celltype)!='rm')],to = celltype_umap)
spe$celltype_umap<-factor(spe$celltype_umap,levels = celltype_umap[!grepl('rm',celltype_umap)])

cols_clst <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
cols_clst <- colorRampPalette(cols_clst)(length(levels(spe$celltype)))
names(cols_clst) <- levels(spe$celltype)

cols_clst_umap <- ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value
cols_clst_umap <- colorRampPalette(cols_clst)(length(levels(spe$celltype_umap)))
names(cols_clst_umap) <- levels(spe$celltype_umap)
spe$pg_cluster_umap<-str_split_fixed(spe$celltype_umap,pattern = ':',n = 2)[,1]
spe$pg_cluster_umap<-factor(spe$pg_cluster_umap %>% as.numeric())

temp <- data.frame(spe$pg_cluster_umap, reducedDim(spe, "UMAP_harmony"))
colnames(temp) <- c('new.id','UMAP_1','UMAP_2')
lc.cent = temp %>% group_by(new.id) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarise_all(median)

hmdata0 <- cbind(t(assay(spe, "exprs"))[,rowData(spe)$use_channel], Anno = spe$celltype) %>%
  as_tibble() %>% group_by(Anno) %>%
  summarise_all(mean) %>% column_to_rownames(var = 'Anno')

cols <- colorRampPalette(rev(brewer.pal(7, 'RdYlBu')))(101)
# Heatmap
temp <- spe$celltype %>% table() %>% as.data.frame()
row_ha = rowAnnotation(counts = anno_barplot(temp$Freq, gp = gpar(fill = cols_clst)))
rownames(hmdata0) <- spe$celltype %>% levels()

dir_create(paste0(output.dir1,'/2.2_Heatmap'))

hmdata0 %>% write.csv(paste(output.dir1, '2.2_Heatmap/heatmap.median.anno.csv', sep = '/'))

### 2.3 marker median heatmap ####
hmdata1 <- cbind(t(assay(spe, "exprs")), Anno = spe$celltype) %>%
  as_tibble() %>% group_by(Anno) %>%
  summarise_all(mean) %>% column_to_rownames(var = 'Anno')

cols <- colorRampPalette(rev(brewer.pal(7, 'RdYlBu')))(101)
# Heatmap
temp <- spe$celltype %>% table() %>% as.data.frame()
row_ha = rowAnnotation(counts = anno_barplot(temp$Freq, gp = gpar(fill = cols_clst)))
rownames(hmdata1) <- spe$celltype %>% levels()

hmdata1 %>% write.csv(paste(output.dir1, '2.2_Heatmap/heatmap.median.anno.all.csv', sep = '/'))


saveRDS(spe, paste(output.dir1, "2_spe.rds", sep = '/'))

outdir11 <- paste0(output.dir1,"/2.8_cluster.marker")
if(! dir.exists(outdir11)) {dir.create(outdir11)}

clu <- unique(spe$celltype) %>% as.character()
group <- as.data.frame(colData(spe))
group <- group[c("roi_id",grep("group",colnames(group),value = T))] %>% unique()


# plot marker expression
if (file.exists(FunctionMarkerCsvDir)) {
  marker <- read.csv(FunctionMarkerCsvDir) %>%
  filter(lineage==1) %>% pull(Antibody )
}else{
  marker <- rownames(assay(spe, "exprs"))[rowData(spe)$use_channel]
}
for (g in grep("group",colnames(group),value = T)) {
  outdir111 <- paste0(outdir11,"/",g)
  if(! dir.exists(outdir111)) {dir.create(outdir111)}
  data <- data.frame(matrix(ncol = 5,nrow = 0))
  colnames(data) <- c("roi_id","group","marker","expresion")
  group1 <- group[c("roi_id",g)]
  colnames(group1)[2] <- "group"
  for (i in clu) {
    spesub <- spe[marker,spe$celltype==i]
    ttest.data <- t(assay(spesub, "exprs")) %>% as.tibble() %>% mutate(roi_id = spesub$roi_id) %>%
      group_by(roi_id) %>% dplyr::summarise_all(median) %>% left_join(group1) %>%
      pivot_longer(cols = !c("roi_id","group"),names_to = "marker",values_to = "expresion")
    ttest.data <- ttest.data[ttest.data$group != '',]

    ttest.data$cluster <- unique(spesub$celltype)
    ttest.data$group <- factor(ttest.data$group,levels = grouporder[g][!grouporder[g]==""])
    ttest.data$marker <- factor(ttest.data$marker,levels = marker)
    data <- rbind(data,ttest.data)

  write.csv(data,paste(outdir111,"expresion.csv",sep = "/"))
}
}


outdir11 <- paste0(output.dir1,'/2.4_echo_roi')
if(!dir.exists(outdir11)) {dir_create(outdir11)}

outdir22 <- paste0(output.dir1,'/2.5_roi.anno')
if(!dir.exists(outdir22)) {dir_create(outdir22)}

## freq table
colData(spe) %>% as_tibble() %>% select(roi_location, celltype) %>% group_by(roi_location) %>% table() %>%
  write.csv(paste0(outdir11,'/freq.csv'))

## roi anno
temp <- cbind(spe@colData,spe@int_colData$spatialCoords) %>% as_tibble() %>% select(roi_location, CellID, celltype,Pos_X,Pos_Y)
colnames(temp)[3] <- 'CellType'
temp.list <- split(temp, temp$roi_location)


lapply(temp.list, function(x) {
  x %>% select(CellID, CellType,Pos_X,Pos_Y)  %>%
    arrange(CellID) %>%
    write.csv(paste(outdir22, paste0(unique(x$roi_location),'.csv'), sep = '/'), row.names = F)
})




marker.df <- t(assay(spe, "exprs")) %>% as.tibble() %>% mutate(roi_location = spe$roi_location,CellID=spe$CellID,Pos_X=spe@int_colData$spatialCoords[,1],Pos_Y=spe@int_colData$spatialCoords[,2])

MarkerName <- colnames(marker.df)
MarkerName <- MarkerName[!MarkerName %in% c("roi_location","CellID","Pos_X","Pos_Y")]

dir_create(paste(output.dir1,'2.6_MarkerNames',sep = '/'))
MarkerName %>% write.csv(paste(output.dir1,'2.6_MarkerNames','MarkerNames.csv',sep = '/'),row.names = F)

roiList <- marker.df$roi_location %>% unique()
dir_create(paste(output.dir1,'2.7_RoiMarker',sep = '/'))
for (roi in roiList){
  marker.df.sub <- marker.df[marker.df$roi_location %in% roi,]
  marker.df.sub %>% write.csv(paste(output.dir1,'2.7_RoiMarker',glue('{roi}.csv'),sep = '/'),row.names = F)

}
print('Finish!')
