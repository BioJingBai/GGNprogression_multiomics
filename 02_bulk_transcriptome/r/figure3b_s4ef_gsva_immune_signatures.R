# -----------------------------------------------------------------------------
# GSVA of the immune-signature gene sets
#
# Purpose   : Compute GSVA scores for the 29 immune signatures in tumour samples and compare the three progression groups.
#
# Inputs    : clinical/ClinicalData.xlsx; bulk_rnaseq/raw_data/mRNA.txt; bulk_rnaseq/reference_data/gsva/29_Immune_signatures.xlsx
#
# Outputs   : figures/GSVA_score_immune29.xlsx (score matrix); figures/GSVA_score_TIL+APC.pdf (group comparison)
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript figure3b_s4ef_gsva_immune_signatures.R
# -----------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

clin<-read.xlsx(file.path(PROJECT_ROOT, "clinical/ClinicalData.xlsx"),sheet = 1)
info<-clin[,c("RNA_Tumor","RNA_Normal","group")]
info_tumor<-info[,c(1,3)]
info_normal<-info[,c(2,3)]
info_tumor<-dplyr::filter(info_tumor,info_tumor$RNA_Tumor!="NG"&info_tumor$group!="NG")
info_normal<-dplyr::filter(info_normal,info_normal$RNA_Normal!="NG"&info_normal$group!="NG")
rna_tumor <- info_tumor$RNA_Tumor[!grepl("NG", info_tumor$RNA_Tumor)]
rna_normal <- info_normal$RNA_Normal[!grepl("NG", info_normal$RNA_Normal)]
mrna_count<-read.table(file.path(PROJECT_ROOT, "bulk_rnaseq/raw_data/mRNA.txt"),header = TRUE)
rownames(mrna_count)<-mrna_count[,1]
mrna_count<-mrna_count[,-1]
mrna_count_tumor <- mrna_count[, rna_tumor]
mrna_count_normal<-mrna_count[, rna_normal]

immune29<-read_excel(file.path(PROJECT_ROOT, "bulk_rnaseq/reference_data/gsva/29_Immune_signatures.xlsx"),col_names = TRUE)
immune29<-immune29[c(3,26)]
# convert to a list; each column is one gene set
immune29_signatures <- lapply(immune29, function(x) {
  # drop NA and empty values, and drop the first row (the pathway name)
  genes <- na.omit(x[-1])
  genes[genes != ""]
})
# use the column names as the list names
names(immune29_signatures) <- immune29[1, ]
# build the object required by GSVA
GSVA_object_immune29<-gsvaParam(as.matrix(mrna_count_tumor),immune29_signatures)
GSVA_score_immune29<-gsva(GSVA_object_immune29)
GSVA_score_immune29<-as.data.frame(GSVA_score_immune29)
writexl::write_xlsx(GSVA_score_immune29,file.path(PROJECT_ROOT, "figures/GSVA_score_immune29.xlsx"))
GSVA_score_immune29<-as.data.frame(t(GSVA_score_immune29))
GSVA_score_immune29$RNA_Tumor<-rownames(GSVA_score_immune29)
GSVA_score_immune29<-left_join(GSVA_score_immune29,info_tumor,by="RNA_Tumor")
# continuous variables
vars <- colnames(GSVA_score_immune29)[!grepl("Pattern|group|T_N|P-value|Correlation|RMSE|RNA_Tumor", colnames(GSVA_score_immune29))]
groups <- c("group")
comparisons <- list(c("S","P"),c("S","MP"),c("P","MP"))
mycolour<-c("S"="#00AFBB","P"="#E7B800","MP"="#FC4E07")
plot_list <- list()
for (i in groups){
  for (var in vars) 
  {
    data2 = GSVA_score_immune29[, c(i, var)]# ensure i and var are treated as character names, not as variables
    data2 <- na.omit(data2)
    colnames(data2)[1] <- "group"
    colnames(data2)[2] <- "value"
    data2$group <- factor(data2$group, levels = c("S", "P", "MP"))
    p <- ggboxplot(data2, x = "group", y = "value", color  ="group", palette = mycolour,add = "jitter") + 
    theme(legend.position = "top") +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_blank())  +
  stat_compare_means(comparisons = comparisons, label = "p.signif",method = "wilcox", p.adjust.method = "BH")+
      stat_compare_means(label.y = 2)
    plot_list[[var]] <- p
  }
  plots <- do.call(grid.arrange, c(plot_list, ncol =2))
  ggsave(file.path(PROJECT_ROOT, "figures/GSVA_score_TIL+APC.pdf"), plots, width = 6, height = 6)
}
