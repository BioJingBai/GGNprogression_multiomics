# -----------------------------------------------------------------------------
# Tumour-microenvironment functional-signature enrichment (DESeq2 + GSEA)
#
# Purpose   : Run DESeq2 between the three GGN progression groups (S / P / MP), and additionally each group vs all other samples, then run GSEA of the ranked gene lists against a curated functional-signature gene set and draw the combined enrichment bar plot.
#
# Inputs    : objects from the upstream R session: cdata (sample metadata), gene_count (gene x sample count matrix), fges (pathway/functional-signature table). No files are read from disk.
#
# Outputs   : per-comparison DESeq2 result CSVs and the combined GSEA bar-plot PDF (gesa_new.pdf) under bulk_rnaseq/diff_expr_stage/data/
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript figure3a_pathway_signature_enrichment.R
# -----------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

setwd(file.path(PROJECT_ROOT, "bulk_rnaseq/diff_expr_stage"))
# Step 1: prepare the sample grouping information
# Assume the original data frame is named original_data (to avoid overwriting it)
sample_info <- cdata |> 
  select(Tumor_ID...6, RNA_Tumor,group) |> 
  filter(RNA_Tumor!= "NG") |> 
  filter(group != "NG") |> 
  filter(group %in% c("S", "P", "MP")) |> 
  distinct()  # ensure there are no duplicate samples
# Step 2: prepare the gene count matrix
# Assume gene_count is the gene expression matrix (rows = genes, columns = samples)
# Make sure the column names of gene_count match sample_info$RNA_Tumor

# Step 3: filter and reorder the count matrix so it matches the sample information
a<-gene_count[[1]]
count_data <- gene_count[, -1]
count_data <- count_data[, colnames(count_data) %in% sample_info$Tumor_ID...6]
count_data <- count_data[, match(sample_info$Tumor_ID...6, colnames(count_data))]
count_data <- as.data.frame(count_data)
rownames(count_data)<-a
# Step 4: create the DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = data.frame(group = sample_info$group, row.names = sample_info$Tumor_ID...6),
  design = ~ group
)

# Step 5: run the differential expression analysis
dds <- DESeq(dds)
# Step 6: define the comparison pairs
comparisons <- list(
  c("S", "P"),
  c("S", "MP"),
  c("P", "MP")
)
# Step 7: run the pairwise comparisons and export the results
for (comp in comparisons) {
  res <- results(
    dds,
    contrast = c("group", comp[1], comp[2]),
    alpha = 0.05
  )
  
  # convert to a data frame and add a gene-name column
  res_df <- as.data.frame(res) |> 
    rownames_to_column("Gene") |> 
    arrange(pvalue)
  # build the file name
  filename <- sprintf("%s_VS_%s.deseq2.csv", comp[1], comp[2])
  # export the results
  write_csv(res_df, file.path(file.path(PROJECT_ROOT, "bulk_rnaseq/diff_expr_stage/data"), filename))
}

# Generate the comparison of each group against all other samples

# Step 5: define the comparison groups (each group vs all other samples)
groups <- c("S", "P", "MP")
# Step 6: run the comparison of each group vs all other samples
for (grp in groups) {
  # create a new grouping variable
  new_group <- ifelse(sample_info$group == grp, grp, "other")
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
  res <- results(
    dds_temp,
    contrast = c("new_group", grp, "other"),
    alpha = 0.05
  )
  # convert to a data frame and add a gene-name column
  res_df <- as.data.frame(res) |> 
    rownames_to_column("Gene") |> 
    arrange(pvalue)
  # build the file name
  filename <- sprintf("%s_vs_other.deseq2.csv", grp)
  # export the results
  write_csv(res_df, file.path(file.path(PROJECT_ROOT, "bulk_rnaseq/diff_expr_stage/data"), filename))
}

## Plotting
fges_long <- reshape::melt(fges,id.vars = "X1",value.name = "gene")
fges_long <- fges_long[!is.na(fges_long$value),c("X1","value")]
colnames(fges_long) <- c("fges","gene")
###tumor------------
df <- data.frame(ID = character(), NES = numeric(), pvalue = numeric(), score = numeric(), comparison = character(), stringsAsFactors = FALSE)
# list all files in the current directory
files <- list.files(file.path(PROJECT_ROOT, "bulk_rnaseq/diff_expr_stage/data"))

# select all tumour differential-expression result files
filtered_files <- files[grep(".*csv$", files)]
# get all comparisons
comparisons <- sapply(strsplit(filtered_files, "\\."), function(x) x[1])
for (i in comparisons) {
  # read the file
  data <- read.csv(paste0(file.path(PROJECT_ROOT, "bulk_rnaseq/diff_expr_stage/data/"),i, ".deseq2.csv"))
  gene <- data$Gene
  gene = bitr(gene, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
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
  # GSEA
  set.seed(123456)
  fges_gsea <- GSEA(geneList, TERM2GENE = fges_long, by = "fgsea",
                    minGSSize = 2, maxGSSize = 1000, pvalueCutoff = 1)
  df <- rbind(df, data.frame(
    ID = fges_gsea$ID,
    NES = fges_gsea$NES,
    pvalue = fges_gsea$pvalue,
    score = abs(log10(fges_gsea$pvalue)) * sign(fges_gsea$NES),
    comparison = i,
    stringsAsFactors = FALSE))
}

# Remove some non-significant results from the tumour comparison
#df <- df[!(df$ID %in% c("Myeloid cells traffic", "Granulocyte traffic","Tumor-associated Macrophages","Th2 signature","Antitumor cytokines","Protumor cytokines","EMT signature","Treg","Treg and Th2 traffic","M1 signature")), ]

# Plot!
# set the position of the red dashed lines
threshold1 <- -log10(0.05)
threshold2 <- log10(0.05)

# create an empty list to store the plot of each facet
plot_list <- list()
group_colors <- c("None" = "gray80", "Others" = "#808285", "S" = "#7570B3", "P" = "steelblue","MP" = "#D96226") 

# loop over every comparison
for (j in comparisons) {
  # select the rows of df for the current comparison
  comp_df <- df[df$comparison == j, ]
  
  # reorder the ID variable by score
  comp_df$ID <- reorder(comp_df$ID, comp_df$score)
  
  comp_df$group <- ifelse(comp_df$score > log10(0.05) & comp_df$score < -log10(0.05), 
                          "None",
                          ifelse(comp_df$score > -log10(0.05),
                                 sub("_.*", "", comp_df$comparison),
                                 sub(".*_", "", comp_df$comparison)))
  
  # create the ggplot object and draw the bar chart
  p <- ggplot(comp_df, aes(x = score, y = ID, fill = group)) +
    geom_bar(stat = "identity",color = "black") +
    theme_bw() +
    scale_fill_manual(values = group_colors) +
    theme(panel.spacing = unit(0.5, "lines"),
          panel.border = element_rect(linewidth = 2),
          panel.grid = element_blank(),
          legend.position = "") +
    labs(title = j,x = "-log10(P-value)", y = "")
  
  # add the red dashed line for the positive threshold
  p <- p + geom_vline(xintercept = threshold1, linetype = "dashed", color = "darkred")
  
  # add the red dashed line for the negative threshold
  p <- p + geom_vline(xintercept = threshold2, linetype = "dashed", color = "darkred")
  
  # add the current facet's plot to the list
  plot_list[[j]] <- p
}

# arrange all facet plots in the specified order
grid_arrange <- function(plot_list) {
  row1 <- c(5,NA,NA,NA)
  row2 <- c(6,3,NA,NA)
  row3 <- c(4,2,1,NA)
  grid.arrange(grobs = plot_list, layout_matrix = rbind(row1,row2,row3), ncol = 3)
}

# combine all facet plots together
plots <- grid_arrange(plot_list)
# export the figure
ggsave(file.path(PROJECT_ROOT, "bulk_rnaseq/diff_expr_stage/data/gesa_new.pdf"), plots, width = 16, height = 15)
