# -----------------------------------------------------------------------------
# TCR clonality across the three growth patterns
#
# Purpose   : Compare tumour TCR clonality (T_clonality) across the three growth-pattern groups (S / P / MP) with a Benjamini-Hochberg-adjusted two-sided Wilcoxon test.
#
# Inputs    : clinical/ClinicalData.xlsx (columns: group, T_clonality)
#
# Outputs   : figures/tcr_clonality.pdf
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript figure2h_tcr_clonality.R
# -----------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

tcr_clonal<-read_excel(file.path(PROJECT_ROOT, "clinical/ClinicalData.xlsx"),sheet = 1)
tcr_clonal<-tcr_clonal[,c("group","T_clonality")]
tcr_clonal<-dplyr::filter(tcr_clonal,tcr_clonal$T_clonality!="NG"&tcr_clonal$T_clonality!="NA"&tcr_clonal$group!="NG")
tcr_clonal$T_clonality<-as.numeric(tcr_clonal$T_clonality)
# define the helper function
#get_adj_p <- function(data, .col, .grp = "Sample", comparisons = NULL,method = "wilcox.test", p.adjust.method = "fdr", p.digits = 3L,symnum.args = NULL, ...) {
#   # Compute p-values
#   comparison.formula <- paste0(.col, "~", .grp) %>%
#     as.formula()
#   pvalues <- ggpubr::compare_means(
#     formula = comparison.formula, data = data,
#     method = method,
#     p.adjust.method = p.adjust.method,
#     ...
#   )
#   
#   # If a comparison list is provided, extract the comparisons of interest for plotting
#   if (!is.null(comparisons)) {
#     pvalues <- purrr::map_df(comparisons, ~ pvalues %>% dplyr::filter(group1 == .x[1] & group2 == .x[2]))
#   }
#   
#   # P-value y coordinates
#   y.max <- data %>%
#     dplyr::pull(.col) %>%
#     max(na.rm = TRUE)
#   p.value.y.coord <- rep(y.max, nrow(pvalues))
#   
#   step.increase <- (1:nrow(pvalues)) * (y.max / 10)
#   p.value.y.coord <- p.value.y.coord + step.increase
#   if (is.null(symnum.args)){
#     symnum.args <- list(cutpoints = c(0, 1e-04, 0.001, 0.01, 
#                                       0.05, 1), symbols = c("****", "***", "**", "*", 
#                                                             "ns"))
#   } 
#   
#   symnum.args$x <- as.numeric(pvalues$p.adj)
#   p.adj.signif <- do.call(stats::symnum, symnum.args) %>% 
#     as.character()
#   pvalues$p.adj.signif = p.adj.signif
#   pvalues <- pvalues %>%
#     dplyr::mutate(
#       y.position = p.value.y.coord,
#       p.adj = format.pval(.data$p.adj, digits = p.digits)
#       
#     )
#   
#   pvalues
# }
#p_adj <- get_adj_p(tcr_clonal,.col = "T_clonality", .grp = "group", p.adjust.method = "BH",comparisons = list(c("S","P"),c("S","MP"),c("P","MP")))
mycolour<-c("S"="#F0D490","P"="#6fA9C8","MP"="#B54030")
tcr_clonal$group <- factor(tcr_clonal$group, levels = c("S", "P", "MP"))
comparisons <- list(c("S","P"),c("S","MP"),c("P","MP"))
p <- ggboxplot(tcr_clonal, x = "group", y = "T_clonality", fill ="group", palette = mycolour,add = "jitter",add.params = list(
        color = "black",   # point colour
        alpha = 0.6,       # transparency
        size = 1.5         # point size
    )) + 
    theme(legend.position = "top") +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x = element_blank())  +
  stat_compare_means(
  label = "p.signif",  # show the adjusted significance labels
  comparisons = comparisons,
  method = "wilcox.test",  # note: use wilcox.test, not wilcox
  p.adjust.method = "BH",  # Benjamini-Hochberg correction
)
pdf(file = file.path(PROJECT_ROOT, "figures/tcr_clonality.pdf"), width =4, height = 5)
print(p)
dev.off()
