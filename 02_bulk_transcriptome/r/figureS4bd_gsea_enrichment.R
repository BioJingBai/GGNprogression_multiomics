# -----------------------------------------------------------------------------
# GSEA enrichment (GO-BP / KEGG / Hallmark)
#
# Purpose   : Run GSEA (GO-BP via gseGO, Hallmark via the MSigDB GMT, KEGG via gseKEGG) on each tumour DEG table and collect the enrichment statistics into three result tables.
#
# Inputs    : bulk_rnaseq/results/Deseq2/DEG_*vsother_tumor.csv (written by script 05); bulk_rnaseq/reference_data/gsea/h.all.v2025.1.Hs.symbols.gmt; org.Hs.eg.db annotation package.
#
# Outputs   : figures/GO_BP_plot.csv, figures/hallmark_plot.csv, figures/KEGG_plot.csv
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript figureS4bd_gsea_enrichment.R
# -----------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

#GO
files <- list.files(file.path(PROJECT_ROOT, "bulk_rnaseq/results/Deseq2"))
comparisons <- c("Svsother","Pvsother","MPvsother")
GO_BP_plot <- data.frame()  # initialise the result data frame
for (i in comparisons) {
  data <- read.csv(paste0(file.path(PROJECT_ROOT, "bulk_rnaseq/results/Deseq2/DEG_"),i, "_tumor.csv"))
  gene <- data$Gene
  gene = bitr(gene, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")# convert using the raw (ENSEMBL) identifiers first
  # de-duplicate genes
  gene <- distinct(gene, SYMBOL, .keep_all = TRUE)
  colnames(gene) <- c("Gene", "SYMBOL")
  # merge data and SYMBOL
  data_all <- data %>% 
    inner_join(gene, by = "Gene")
  # sort genes by logFC and p-value from high to low
  data_all_sort <- data_all %>% 
    arrange(desc(log2FoldChange))
  # extract the fold changes from largest to smallest
  geneList <- data_all_sort$log2FoldChange
  # attach the matching SYMBOL to each extracted fold change
  names(geneList) <- data_all_sort$SYMBOL
  # make sure geneList contains no NA, otherwise GSEA errors out
  geneList <- na.omit(geneList)
  set.seed(123456)
  GO_BP_kk_entrez <- gseGO(geneList  = geneList,
                      ont          = "BP",  # "BP", "MF", "CC" or "ALL"
                      OrgDb        = "org.Hs.eg.db",# human: org.Hs.eg.db; mouse: org.Mm.eg.db
                      keyType      = "SYMBOL",
                      pvalueCutoff = 0.25)   # this is actually the padj threshold and can be adjusted
  GO_BP_kk <- DOSE::setReadable(GO_BP_kk_entrez, 
                           OrgDb="org.Hs.eg.db",
                           keyType="SYMBOL")# convert IDs 
  GO_BP_plot <- rbind(GO_BP_plot, data.frame(
    ID = GO_BP_kk$Description,
    NES = GO_BP_kk$NES,
    pvalue = GO_BP_kk$pvalue,
    score = abs(log10(GO_BP_kk$pvalue)) * sign(GO_BP_kk$NES),
    comparison = i,
    stringsAsFactors = FALSE))
}
write_csv(GO_BP_plot,file.path(PROJECT_ROOT, "figures/GO_BP_plot.csv"))

#Hallmark
files <- list.files(file.path(PROJECT_ROOT, "bulk_rnaseq/results/Deseq2"))
h<-read.gmt(file.path(PROJECT_ROOT, "bulk_rnaseq/reference_data/gsea/h.all.v2025.1.Hs.symbols.gmt"))
comparisons <- c("Svsother","Pvsother","MPvsother")
h_plot <- data.frame()  # initialise the result data frame
for (i in comparisons) {
  data <- read.csv(paste0(file.path(PROJECT_ROOT, "bulk_rnaseq/results/Deseq2/DEG_"),i, "_tumor.csv"))
  gene <- data$Gene
  gene = bitr(gene, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")# convert using the raw (ENSEMBL) identifiers first
  # de-duplicate genes
  gene <- distinct(gene, SYMBOL, .keep_all = TRUE)
  colnames(gene) <- c("Gene", "SYMBOL")
  # merge data and SYMBOL
  data_all <- data %>% 
    inner_join(gene, by = "Gene")
  # sort genes by logFC and p-value from high to low
  data_all_sort <- data_all %>% 
    arrange(desc(log2FoldChange))
  # extract the fold changes from largest to smallest
  geneList <- data_all_sort$log2FoldChange
  # attach the matching SYMBOL to each extracted fold change
  names(geneList) <- data_all_sort$SYMBOL
  # make sure geneList contains no NA, otherwise GSEA errors out
  geneList <- na.omit(geneList)
  set.seed(123456)
  fges_gsea <- GSEA(geneList, TERM2GENE = h, by = "fgsea",
                    minGSSize = 2, maxGSSize = 1000, pvalueCutoff = 1)
  fges_gsea<-fges_gsea[fges_gsea$pvalue < 0.05, ]# keep only significant results before merging into df
  h_plot <- rbind(h_plot, data.frame(
    ID = fges_gsea$ID,
    NES = fges_gsea$NES,
    pvalue = fges_gsea$pvalue,
    score = abs(log10(fges_gsea$pvalue)) * sign(fges_gsea$NES),
    comparison = i,
    stringsAsFactors = FALSE))
}
write_csv(h_plot,file.path(PROJECT_ROOT, "figures/hallmark_plot.csv"))

#KEGG
options(timeout = 300)
files <- list.files(file.path(PROJECT_ROOT, "bulk_rnaseq/results/Deseq2"))
comparisons <- c("Svsother","Pvsother","MPvsother")
KEGG_plot <- data.frame()  # initialise the result data frame
for (i in comparisons) {
  data <- read.csv(paste0(file.path(PROJECT_ROOT, "bulk_rnaseq/results/Deseq2/DEG_"),i, "_tumor.csv"))
   gene <- data$Gene
  gene = bitr(gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")# convert using the raw (ENSEMBL) identifiers first
  # de-duplicate genes
  gene <- distinct(gene, ENTREZID, .keep_all = TRUE)
  colnames(gene) <- c("Gene", "ENTREZID")
  # merge data and SYMBOL
  data_all <- data %>% 
    inner_join(gene, by = "Gene")
  # sort genes by logFC and p-value from high to low
  data_all_sort <- data_all %>% 
    arrange(desc(log2FoldChange))
  # extract the fold changes from largest to smallest
  geneList <- data_all_sort$log2FoldChange
  # attach the matching SYMBOL to each extracted fold change
  names(geneList) <- data_all_sort$ENTREZID
  # make sure geneList contains no NA, otherwise GSEA errors out
  geneList <- na.omit(geneList)
  set.seed(123456)
  KEGG_kk_entrez <- gseKEGG(geneList   = geneList,
                          organism     = "hsa", # human: hsa; mouse: mmu
                          pvalueCutoff = 0.25)  # this is actually the padj threshold and can be adjusted 
  KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                             OrgDb="org.Hs.eg.db",
                             keyType='ENTREZID')# convert IDs            
  KEGG_plot <- rbind(KEGG_plot, data.frame(
    ID = KEGG_kk$Description,
    NES = KEGG_kk$NES,
    pvalue = KEGG_kk$pvalue,
    score = abs(log10(KEGG_kk$pvalue)) * sign(KEGG_kk$NES),
    comparison = i,
    stringsAsFactors = FALSE))
}
write_csv(KEGG_plot,file.path(PROJECT_ROOT, "figures/KEGG_plot.csv"))
