# -----------------------------------------------------------------------------
# DESeq2 differential expression in tumour samples
#
# Purpose   : Run DESeq2 for the three progression groups (S / P / MP) and additionally each group vs all other samples, writing the differential-expression tables consumed by the GSEA scripts.
#
# Inputs    : objects from the upstream R session: info (sample metadata, created by script 03) and rna_count_tumor (tumour count matrix). No files are read from disk.
#
# Outputs   : DEG_SvsP_tumor.csv, DEG_SvsMP_tumor.csv, DEG_PvsMP_tumor.csv and DEG_<group>vsother_tumor.csv under bulk_rnaseq/results/Deseq2/
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript figureS4bd_deseq2_tumor.R
# -----------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

# Tumor
group<- factor(c(rep("S",8),rep("P",15),rep("MP",19)))
coldata_tumor<-info[,c(1,3)]
coldata_tumor<-dplyr::filter(coldata_tumor,coldata_tumor$RNA_Tumor!="NG"&coldata_tumor$group!="NG")
rownames(coldata_tumor) <- colnames(rna_count_tumor)
dds <- DESeqDataSetFromMatrix(countData = rna_count_tumor,
                              colData = coldata_tumor,
                              design = ~ group)
head(dds) # inspect the constructed matrix
dds <- DESeq(dds)
# define the comparisons
res_SvsP <- results(dds, contrast = c("group", "S", "P"))
res_SvsMP <- results(dds, contrast = c("group", "S", "MP"))
res_PvsMP <- results(dds, contrast = c("group", "P", "MP"))
# convert to a data.frame and sort by p-value
res_SvsP_df  <- as.data.frame(res_SvsP[order(res_SvsP$pvalue), ]) |>
  rownames_to_column("Gene")
res_SvsMP_df <- as.data.frame(res_SvsMP[order(res_SvsMP$pvalue), ])|>
  rownames_to_column("Gene")
res_PvsMP_df <- as.data.frame(res_PvsMP[order(res_PvsMP$pvalue), ])|>
  rownames_to_column("Gene")
# export CSV files
write.csv(res_SvsP_df,  file = file.path(PROJECT_ROOT, "bulk_rnaseq/results/Deseq2/DEG_SvsP_tumor.csv"))
write.csv(res_SvsMP_df, file = file.path(PROJECT_ROOT, "bulk_rnaseq/results/Deseq2/DEG_SvsMP_tumor.csv"))
write.csv(res_PvsMP_df, file = file.path(PROJECT_ROOT, "bulk_rnaseq/results/Deseq2/DEG_PvsMP_tumor.csv"))

# each group vs the other two groups
groups <- c("S", "P", "MP")
# Step 6: run the comparison of each group vs all other samples
for (grp in groups) {
  # create a new grouping variable
  new_group <- ifelse(coldata_tumor$group == grp, grp, "other")
  # create a temporary dds object
  dds_temp <- dds
  colData(dds_temp)$new_group <- factor(new_group)
  # set the reference level to "other"
  dds_temp$new_group <- relevel(dds_temp$new_group, ref = "other")
  # update the design formula
  design(dds_temp) <- ~ new_group
  # re-run DESeq (reusing the existing size factors)
  dds_temp <- DESeq(dds_temp, 
                   test = "Wald",
                   fitType = "parametric",
                   sfType = "ratio",
                   betaPrior = TRUE,
                   modelMatrixType = "expanded",
                   useT = FALSE,
                   minmu = 1e-6,
                   parallel = FALSE)
  # extract the results
  res <- results(dds_temp,contrast = c("new_group", grp, "other"),alpha = 0.05)
  # convert to a data frame and add a gene-name column
  res_df <- as.data.frame(res) |> 
    rownames_to_column("Gene") |> 
    arrange(pvalue)
  # build the file name
  filename <- sprintf("DEG_%svsother_tumor.csv", grp)
  # export the results
  write_csv(res_df, file.path(file.path(PROJECT_ROOT, "bulk_rnaseq/results/Deseq2"), filename))
}
