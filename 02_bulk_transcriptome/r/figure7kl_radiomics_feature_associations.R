# -----------------------------------------------------------------------------
# Association of the retained radiomics features with immune-related genes
#
# Purpose   : Build the bulk-RNA feature table - GSVA scores of the Hallmark, TLS and FgSE gene sets merged with immune-cell infiltration and ESTIMATE scores - that is associated with the retained radiomics features.
#
# Inputs    : clinical/ClinicalData.xlsx; bulk_rnaseq/raw_data/mRNA.txt; bulk_rnaseq/diff_expr_stage/20250314/fgse.xlsx; bulk_rnaseq/immune_infiltration/20250313/immune_score_normlization.tsv; bulk_rnaseq/immune_infiltration/20250313/estimate_results.csv
#
# Outputs   : bulk_rnaseq/features_bulk.csv (full feature table); bulk_rnaseq/features_bulk_scores.csv (score-only table)
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript figure7kl_radiomics_feature_associations.R
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

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

TLS_signatures<-list(cxcl13="CXCL3",chem12=c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13"),Tfh=c("CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1"),Th1_B=c("CD4","CCR5","CXCR3","CSF2","IGSF6","IL-2RA","CD38","CD40","CD5","MS4A1","SDC1","GFI1","IL-1R1","IL-1R2","IL-10","CCL20","TRAF6","STAT5A"),PlasmaCell=c("TNFRSF17","IGJ"),TLS_imprint=c("IGHA1","IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGHM","IGKC","IGLC1","IGLC2","IGLC3","JCHAIN","CD79A","FCRL5","MZB1","SSR4","XBP1","TRBC2","IL-7R","CXCL12","LUM","C1QA","C7","CD52","APOE","PTLP","PTGDS","PIM2","DERL3"),TLS_DSP=c("CYTIP","CXCL9","CD79A","LIMD2","LSP1","SELL","BIRC3","PLEK","TRBC1","POU2AF1","RCSD1","ARHGDIB","LTB","TAGAP","CD52","WIPF1","CORO1A","PTPRC","LCP1","PTPRCAP"),TLS_genes=c("CXCL3","CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13","CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1","CD4","CCR5","CXCR3","CSF2","IGSF6","IL-2RA","CD38","CD40","CD5","MS4A1","SDC1","GFI1","IL-1R1","IL-1R2","IL-10","CCL20","TRAF6","STAT5A","TNFRSF17","IGJ","IGHA1","IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGHM","IGKC","IGLC1","IGLC2","IGLC3","JCHAIN","CD79A","FCRL5","MZB1","SSR4","XBP1","TRBC2","IL-7R","CXCL12","LUM","C1QA","C7","CD52","APOE","PTLP","PTGDS","PIM2","DERL3","CYTIP","CXCL9","CD79A","LIMD2","LSP1","SELL","BIRC3","PLEK","TRBC1","POU2AF1","RCSD1","ARHGDIB","LTB","TAGAP","CD52","WIPF1","CORO1A","PTPRC","LCP1","PTPRCAP"))

fgse<-read.xlsx(file.path(PROJECT_ROOT, "bulk_rnaseq/diff_expr_stage/20250314/fgse.xlsx"),colNames = FALSE)
colnames(fgse)<-c("pathway","gene")
fgse_list <- split(fgse$gene, fgse$pathway)
gene_sets <- c(
  hallmark_list,
  TLS_signatures,
  fgse_list
)
gsva_param <- gsvaParam(
  expr = as.matrix(mrna_count_tumor),
  geneSets = gene_sets,
  kcdf = "Gaussian",   # bulk RNA must use Gaussian
  maxDiff = TRUE
)
gsva_res <- gsva(gsva_param)
gsva_df <- as.data.frame(t(gsva_res))
gsva_df$ID <- rownames(gsva_df)

mrna<-as.data.frame(t(mrna_count_tumor))
mrna$ID <- rownames(mrna)
immunecell<-fread(file.path(PROJECT_ROOT, "bulk_rnaseq/immune_infiltration/20250313/immune_score_normlization.tsv"))
immunescore<-read_csv(file.path(PROJECT_ROOT, "bulk_rnaseq/immune_infiltration/20250313/estimate_results.csv"))
immunecell<-as.data.frame(t(immunecell))
colnames(immunecell)<-immunecell[1,]
immunecell<-immunecell[-1,]
immunecell$ID <- paste0("X", rownames(immunecell))
final_df <- immunescore %>%
  dplyr::left_join(immunecell, by = "ID") %>%
  dplyr::left_join(gsva_df, by = "ID") %>%
  dplyr::left_join(
   mrna,by = "ID")
final_df1 <- immunescore %>%
  dplyr::left_join(immunecell, by = "ID") %>%
  dplyr::left_join(gsva_df, by = "ID") 
write_csv(final_df,file.path(PROJECT_ROOT, "bulk_rnaseq/features_bulk.csv"))
write_csv(final_df1,file.path(PROJECT_ROOT, "bulk_rnaseq/features_bulk_scores.csv"))
