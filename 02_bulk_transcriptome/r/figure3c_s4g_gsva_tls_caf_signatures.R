# -----------------------------------------------------------------------------
# GSVA of the TLS / CAF signature gene sets
#
# Purpose   : Compute GSVA scores for the CAF/TLS gene sets in tumour samples and compare the three progression groups.
#
# Inputs    : bulk_rnaseq/reference_data/gsva/CAF_and_TLS.xlsx; objects from the upstream R session: mrna_count_tumor and info_tumor (created by script 03).
#
# Outputs   : figures/GSVA_score_caf_tls.xlsx (score matrix); figures/GSVA_score_TLS.pdf (group comparison)
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript figure3c_s4g_gsva_tls_caf_signatures.R
# -----------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

caf_tls<-read_excel(file.path(PROJECT_ROOT, "bulk_rnaseq/reference_data/gsva/CAF_and_TLS.xlsx"),col_names = FALSE)
caf_tls<-caf_tls[c(3:6),]
caf_tls_signatures <- setNames(
  lapply(1:nrow(caf_tls), function(i) {
    row_genes <- as.character(unlist(caf_tls[i, -1]))
    row_genes[!is.na(row_genes) & row_genes != ""]
  }),
  caf_tls[[1]]
)
# build the object required by GSVA
GSVA_object_caf_tls<-gsvaParam(as.matrix(mrna_count_tumor),caf_tls_signatures)
GSVA_score_caf_tls<-gsva(GSVA_object_caf_tls)
GSVA_score_caf_tls<-as.data.frame(GSVA_score_caf_tls)
writexl::write_xlsx(GSVA_score_caf_tls,file.path(PROJECT_ROOT, "figures/GSVA_score_caf_tls.xlsx"))

comparisons <- list( c("S", "P"), c("P", "MP"), c("S", "MP") )

#tls
GSVA_score_TLS<-GSVA_score_caf_tls
GSVA_score_TLS<-as.data.frame(t(GSVA_score_TLS))
GSVA_score_TLS$RNA_Tumor<-rownames(GSVA_score_TLS)
GSVA_score_TLS<-left_join(GSVA_score_TLS,info_tumor,by="RNA_Tumor")
# continuous variables
vars <- colnames(GSVA_score_TLS)[!grepl("Pattern|group|T_N|P-value|Correlation|RMSE|RNA_Tumor", colnames(GSVA_score_TLS))]
groups <- c("group")
comparisons <- list(c("S","P"),c("S","MP"),c("P","MP"))
plot_list <- list()
for (i in groups){
  for (var in vars) 
  {
    data2 = GSVA_score_TLS[, c(i, var)]# ensure i and var are treated as character names, not as variables
    data2 <- na.omit(data2)
    colnames(data2)[1] <- "group"
    colnames(data2)[2] <- "value"
    data2$group <- factor(data2$group, levels = c("S", "P", "MP"))
    p <- ggviolin(data2, x = "group", y = "value", fill = "group",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white")) +
      labs(x = "", y =  var) +
      theme(axis.text=element_text(face="bold", size=14,color="black"),
            axis.title.y=element_text(face="bold", size=14,color ="black"),
            axis.text.x = element_text(),
            strip.text = element_text(size = 14,face = "bold"),
            strip.background = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none" ) +
   stat_compare_means(comparisons = comparisons, label = "p.signif",method = "wilcox", p.adjust.method = "BH")+
      stat_compare_means(label.y = 2)
    plot_list[[var]] <- p
  }
  plots <- do.call(grid.arrange, c(plot_list, ncol =2))
  ggsave(file.path(PROJECT_ROOT, "figures/GSVA_score_TLS.pdf"), plots, width = 8, height = 8)
}
