# -----------------------------------------------------------------------------
# Differential cell-cell interactions (imcRtools classic / histocat / patch)
#
# Purpose   : Build the k-NN spatial graph and test for significant
#             interaction / avoidance between every pair of cell types with
#             imcRtools::testInteractions / countInteractions, then write the
#             interaction count matrices, the per-group fold-change matrices and
#             the significance tables.
#
# Inputs    : <PROJECT_ROOT>/imc/out/1_Cluster/1_spe.rds
#
# Outputs   : <PROJECT_ROOT>/imc/out/5_Interaction/spe_Interaction.rds
#             <PROJECT_ROOT>/imc/out/5_Interaction/spe_Interaction_out.rds
#             <PROJECT_ROOT>/imc/out/5_Interaction/spe_Interaction_ct.rds
#             <PROJECT_ROOT>/imc/out/5_Interaction/All/*.csv
#             <PROJECT_ROOT>/imc/out/5_Interaction/<group>/*.csv
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 05_cell_cell_interactions.R
# -----------------------------------------------------------------------------

library(RColorBrewer)
library(wesanderson)
library(paletteer)
library(ggplot2)
library(ComplexHeatmap)
library(ggthemes)
library(ggpubr)
library(plyr)
library(ggbeeswarm)
library(argparser)
#library(scales)
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
library(glue)

set.seed(100)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

rds <- file.path(PROJECT_ROOT, "imc", "out", "1_Cluster", "1_spe.rds")
spe <- readRDS(rds)
method <- 'classic'
graphType <- 'knn'

# output dir
output.dir <- file.path(PROJECT_ROOT, "imc", "out")
if(! dir.exists(output.dir)) {dir.create(output.dir)}

group.stat.list <- list(c("T","N"),c("T","MT"),c("N","MT"))
knn_k = 20

set.seed(220705)

if (graphType=='knn'){
  spe <- buildSpatialGraph(spe, img_id = "roi_id", type = "knn", k = knn_k)
  spe <- aggregateNeighbors(spe, colPairName = "knn_interaction_graph",
                            aggregate_by = "metadata", count_by = "celltype")
  cPName <- "knn_interaction_graph"
}else if (graphType=='expansion'){
  spe <- buildSpatialGraph(spe, img_id = "roi_id", type = "expansion", threshold = knn_k)
  spe <- aggregateNeighbors(spe, colPairName = "expansion_interaction_graph",
                            aggregate_by = "metadata", count_by = "celltype")
  cPName <- "expansion_interaction_graph"
}



## 5. Interaction ############################################################
output.dir5 <-  paste(output.dir, '5_Interaction', sep = '/')
if(! dir.exists(output.dir5)) {dir.create(output.dir5)}
#method<-'classic'
if(method %in% c('classic','histocat')) {
  out <- testInteractions(spe,
                          group_by = "roi_id",
                          label = "celltype",
                          method = method,
                          colPairName = cPName,
                          BPPARAM = SerialParam(RNGseed = 221029)) %>% as_tibble()
  ct <- countInteractions(
    spe,
    group_by = "roi_id",
    label = "celltype",
    method = method,
    colPairName = cPName) %>% as_tibble()
}

if(method == 'patch') {
  out <- testInteractions(spe,
                          group_by = "roi_id",
                          label = "celltype",
                          method = method,
                          patch_size = 3,
                          colPairName = cPName,
                          BPPARAM = SerialParam(RNGseed = 221029)) %>% as_tibble()
  ct <- countInteractions(
    spe,
    group_by = "roi_id",
    label = "celltype",
    method = method,
    patch_size = 3,
    colPairName = cPName) %>% as_tibble()
}

saveRDS(object = spe,file = paste(output.dir5, 'spe_Interaction.rds', sep = '/'))
saveRDS(object = out,file = paste(output.dir5, 'spe_Interaction_out.rds', sep = '/'))
saveRDS(object = ct,file = paste(output.dir5, 'spe_Interaction_ct.rds', sep = '/'))

colnames(ct)[1] <- 'roi_id'
idx <- ct$from_label %>% unique() %>% str_split(., pattern = ':', simplify = T)
idx <- idx[,1] %>% as.integer() %>% order
ct$from_label <- factor(ct$from_label, levels = levels(ct$from_label)[idx])
idx <- ct$to_label %>% unique() %>% str_split(., pattern = ':', simplify = T)
idx <- idx[,1] %>% as.integer() %>% order
ct$to_label <- factor(ct$to_label, levels = levels(ct$to_label)[idx])

meta.info <- colData(spe)[, grepl("group|roi_id", colnames(colData(spe)))] %>% as_tibble() %>% unique()

dir.create(paste0(output.dir5,'/All/'))

ct %>% left_join(meta.info) %>% write.csv(paste(output.dir5, '/All/Interaction_counts.csv', sep = '/'))

cols_hm <- colorRampPalette(brewer.pal(7, 'Greens'))(100)

ct %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  dplyr::summarize(mean_ct = mean(ct, na.rm = TRUE)) %>%
  pivot_wider(names_from = from_label, values_from = mean_ct) %>%
  column_to_rownames(var = 'to_label') %>%
  as.matrix() %>%
  write.csv(paste(output.dir5, paste0('/All/Interaction_counts_mat.csv'), sep = '/'))

temp <- ct %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  dplyr::summarize(mean_ct = mean(ct, na.rm = TRUE))


# meta.info <- as.data.frame(colData(spe)[, grepl("group", colnames(colData(spe)))])
meta.info <- as.data.frame(colData(spe)[grepl("group", colnames(colData(spe)))])
if(ncol(meta.info) == 1) {colnames(meta.info) = 'group1'}
hmdata <- cbind(meta.info, roi_id=spe$roi_id) %>% unique() %>% right_join(as.data.frame(ct))
group.num <- colnames(hmdata)[grep('group', colnames(hmdata))]

for (g in group.num) {
  stat.dir <- paste(output.dir5, g, sep = '/')
  if(! dir.exists(stat.dir)) {dir.create(stat.dir)}

  temp <- hmdata %>% as_tibble() %>%
    select(from_label, to_label,roi_id, g, ct)
  colnames(temp)[4] <- 'group'

  ## t.test
  if (group.stat.list[g][[1]][['stat']] == 't.test') {
    stat.label <- group.stat.list[g][[1]][['stat.label']]
    for(stat.num in 1:length(stat.label)) {
      group.label <- stat.label[[stat.num]]

      hmdata0 <- temp %>% filter(group == group.label[1])
      m0 <- hmdata0 %>% as_tibble() %>%
        group_by(from_label, to_label) %>%
        dplyr::summarize(mean_ct = mean(ct, na.rm = TRUE)) %>%
        pivot_wider(names_from = from_label, values_from = mean_ct) %>%
        column_to_rownames(var = 'to_label') %>%
        as.matrix()
      m0[is.nan(m0)] <- 0

      hmdata1 <- temp %>% filter(group == group.label[2])
      m1 <- hmdata1 %>% as_tibble() %>%
        group_by(from_label, to_label) %>%
        dplyr::summarize(mean_ct = mean(ct, na.rm = TRUE)) %>%
        pivot_wider(names_from = from_label, values_from = mean_ct) %>%
        column_to_rownames(var = 'to_label') %>%
        as.matrix()
      m1[is.nan(m1)] <- 0

      tdata <- temp %>%
        group_by(from_label, to_label,roi_id, group) %>%
        dplyr::summarize(mean_ct = mean(ct, na.rm = TRUE))
      tdata$mean_ct[is.nan(tdata$mean_ct)] <- 0

      p.out <- t(data.frame(1:3))[0,]
      colnames(p.out) <- c('from','to','p.val')
      for (f in unique(tdata$from_label)) {
        for (t in unique(tdata$to_label)) {
          print(c(f,t))
          g0 <- tdata %>% as.data.frame %>% filter(from_label == f, to_label==t, group == group.label[1]) %>% pull(mean_ct)
          g1 <- tdata %>% as.data.frame %>% filter(from_label == f, to_label==t, group == group.label[2]) %>% pull(mean_ct)
          p.val <- t.test(g0,g1)$p.value
          p.out <- rbind(p.out, c(f,t,p.val))
        }
      }
      p.mat <- as.tibble(p.out) %>%
        pivot_wider(names_from = from, values_from = p.val) %>%
        column_to_rownames(var = 'to')

      p.mat <- p.mat %>%
        mutate(across(everything(), as.numeric))
       p.mat[is.na(p.mat)]=1

      m0 %>% write.csv(paste(stat.dir, paste0('Interaction_counts_',g, '_',paste(group.label, collapse = 'vs'),'_mat.', group.label[1],'.csv'), sep = '/'))
      m1 %>% write.csv(paste(stat.dir, paste0('Interaction_counts_',g, '_',paste(group.label, collapse = 'vs'),'_mat.', group.label[2],'.csv'), sep = '/'))
      p.mat %>% write.csv(paste(stat.dir, paste0('Interaction_counts_',g, '_',paste(group.label, collapse = 'vs'),'_mat.','pval','.csv'), sep = '/'))


      p.out <- t(data.frame(1:4))[0,]
      colnames(p.out) <- c('from','to','p.val','log2FC')

      for (f in unique(tdata$from_label)) {
        for (t in unique(tdata$to_label)) {
          print(c(f,t))
          g0 <- tdata %>% as.data.frame %>% filter(from_label == f, to_label==t, group == group.label[1]) %>% pull(mean_ct)
          g1 <- tdata %>% as.data.frame %>% filter(from_label == f, to_label==t, group == group.label[2]) %>% pull(mean_ct)
          mg0 <- mean(g0)
          mg1 <- mean(g1)
          if (mg0==0 & mg1==0){
            log2FC <- log2(1)
          }else if (mg0!=0&mg1==0){
            log2FC <- 'max'
          }else if (mg0==0&mg1!=0){
            log2FC <- '-max'
          }else{
            log2FC <- log2(mg0/mg1)
          }
          p.val <- t.test(g0,g1)$p.value
          p.out <- rbind(p.out, c(f,t,p.val,log2FC))
        }
      }
      colnames(p.out)[1:2] <- c('from_label','to_label')
      p.out.1 <- p.out %>% as.data.frame()
      minFC <- p.out.1$log2FC[p.out.1$log2FC!='max'& p.out.1$log2FC!='-max'] %>% as.numeric() %>% min
      maxFC <- p.out.1$log2FC[p.out.1$log2FC!='max'& p.out.1$log2FC!='-max'] %>% as.numeric() %>% max

      tdata.sub <- tdata %>% filter(group!='')

      mat <- tdata.sub %>% as_tibble() %>%
        select(from_label, to_label,roi_id, mean_ct)
      mat <- mat  %>%
        group_by(from_label,to_label) %>% dplyr::summarise_at("mean_ct", mean) %>% replace(is.na(.),0)
      mat <- mat %>% left_join(p.out %>% as_tibble())
      mat$log2FC[mat$log2FC == 'max'] <- maxFC
      mat$log2FC[mat$log2FC == '-max'] <- minFC
      mat$log2FC <- as.numeric(mat$log2FC)
      mat$p.val <- as.numeric(mat$p.val)

      mat$label <- ifelse((abs(mat$log2FC)>1 & mat$p.val<0.05),"|log2(FoldChange)|>1\np<0.05","NotSignificant")

      write.csv(mat,file = paste(stat.dir, paste0('5_Interaction_counts_',g, '_',paste(group.label, collapse = 'vs'),'.csv'), sep = '/'))

    }
  }

}


meta.info <- as.data.frame( colData(spe)[grepl("group", colnames(colData(spe)))] )
colnames(ct)[1] <- 'roi_id'
colnames(out)[1] <- 'roi_id'
ct.hmdata <- cbind(meta.info, roi_id=spe$roi_id) %>% unique() %>% right_join(as.data.frame(ct))
out.hmdata <- cbind(meta.info, roi_id=spe$roi_id) %>% unique() %>% right_join(as.data.frame(out))

if (ct.hmdata[is.na(ct.hmdata$ct),]$ct %>% length() != 0){
  ct.hmdata[is.na(ct.hmdata$ct),]$ct <- 0
}
if (out.hmdata[is.na(out.hmdata$ct),]$ct %>% length() != 0){
  out.hmdata[is.na(out.hmdata$ct),]$ct <- 0
}


group.num <- colnames(hmdata)[grep('group', colnames(hmdata))]




####all out&ct######
mat <- ct.hmdata %>% as_tibble() %>%select(from_label, to_label,roi_id, ct)
#colnames(mat)[4] <- 'group'
mat <- mat %>% group_by(from_label,to_label) %>% dplyr::summarise_at("ct", mean) %>% replace(is.na(.),0)
san <- out.hmdata %>% as_tibble() %>%select(from_label, to_label,roi_id, sigval)
#colnames(san)[4] <- 'group'
#san <- san %>% filter(group==group.label)
nt <- length(unique(san$roi_id))
attr <- dplyr::count(san[san$sigval==1,],from_label,to_label) %>% mutate(attr=n/nt)
avod <- dplyr::count(san[san$sigval==-1,],from_label,to_label) %>% mutate(avod=n/nt)
mat <- left_join(mat, attr, by=c("from_label", "to_label")) %>% left_join(., avod, by=c("from_label", "to_label")) %>% replace(is.na(.),0)
mat$group <- ifelse(mat$attr>mat$avod,"Interaction","Avoidance")
mat[mat$attr==mat$avod,"group"] <- "Tie"
mat$percent <- apply(mat[,c("attr","avod")],1,function(x)max(x))
#mat <- mat %>%  mutate(from_label=factor(from_label, levels = rev(order)))  %>%  mutate(to_label=factor(to_label, levels = order))
mat$to_label<-factor(mat$to_label,levels = rev(unique(mat$to_label)))

####heatmap
mat$to_label<-factor(mat$to_label,levels = rev(unique(mat$to_label)))

mat$`Interaction score`<-mat$n.x-mat$n.y
max_color_num<-max(abs(min(mat$`Interaction score`)),max(mat$`Interaction score`))
colors<- brewer.pal(11, "RdBu")

write.csv(mat,file=paste(output.dir5,  paste0('/All/Interaction_out_ct_all.csv'),
                         sep = '/'))





#####split group.  out&ct########
for (g in group.num) {
  stat.dir <- paste(output.dir5, g, sep = '/')
  if(! dir.exists(stat.dir)) {dir.create(stat.dir)}
  stat.label <- group.stat.list[g][[1]][['stat.label']]

  for(stat.num in 1:length(stat.label)) {
    for (group.num in 1:length(stat.label[[stat.num]])){
      group.label <- stat.label[[stat.num]][group.num]

      mat <- ct.hmdata %>% as_tibble() %>%
        select(from_label, to_label,roi_id, g, ct)
      colnames(mat)[4] <- 'group'
      mat <- mat %>% filter(group==group.label) %>%
        group_by(from_label,to_label) %>% dplyr::summarise_at("ct", mean) %>% replace(is.na(.),0)

      san <- out.hmdata %>% as_tibble() %>%
        select(from_label, to_label,roi_id, g, sigval)
      colnames(san)[4] <- 'group'
      san <- san %>% filter(group==group.label)
      nt <- length(unique(san$roi_id))
      attr <- dplyr::count(san[san$sigval==1,],from_label,to_label) %>% mutate(attr=n/nt)
      avod <- dplyr::count(san[san$sigval==-1,],from_label,to_label) %>% mutate(avod=n/nt)
      mat <- left_join(mat, attr, by=c("from_label", "to_label")) %>%
        left_join(., avod, by=c("from_label", "to_label")) %>% replace(is.na(.),0)
      mat$group <- ifelse(mat$attr>mat$avod,"Interaction","Avoidance")
      mat[mat$attr==mat$avod,"group"] <- "Tie"
      mat$percent <- apply(mat[,c("attr","avod")],1,function(x)max(x))
      #mat <- mat %>%  mutate(from_label=factor(from_label, levels = rev(order)))  %>%  mutate(to_label=factor(to_label, levels = order))

      mat$to_label<-factor(mat$to_label,levels = rev(unique(mat$to_label)))

      write.csv(mat,file=paste(stat.dir,  paste0('Interaction_out_ct_',g,'_',group.label,'.csv'), sep = '/'))

      mat$`Interaction score`<-mat$n.x-mat$n.y
      max_color_num<-max(abs(min(mat$`Interaction score`)),max(mat$`Interaction score`))
      colors<- brewer.pal(11, "RdBu")


    }
  }
}



print('Finish!')
