# -----------------------------------------------------------------------------
# TLS gene-set score in the tumour core
#
# Purpose   : Score every CosMx SMI cell for a TLS gene set (AddModuleScore) and
#             compare the score between the S, P and MP groups inside the tumour
#             core (tissue == "T"), with Kruskal-Wallis and pairwise Wilcoxon
#             tests.
#
# Inputs    : <PROJECT_ROOT>/processed/smi_annotated.qs
#             <PROJECT_ROOT>/raw_data/TLS.xlsx   (TLS gene list, row 3, cols 2:30)
#
# Outputs   : <PROJECT_ROOT>/figures/TLS_score_tumor_core.pdf
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 05_figureS15c_tls_score.R
# -----------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Dependencies
## ---------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(readxl)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(qs)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
dir.create(file.path(PROJECT_ROOT, "figures"), recursive = TRUE, showWarnings = FALSE)

## Load the annotated object produced by 00_setup_and_annotation.R
smi <- qread(file.path(PROJECT_ROOT, "processed", "smi_annotated.qs"))


## ---------------------------------------------------------------------------
## TLS gene set and module score
## ---------------------------------------------------------------------------
tls_gene <- read_excel(file.path(PROJECT_ROOT, "raw_data", "TLS.xlsx"))
tls_genes <- as.character(unlist(tls_gene[3, 2:30]))
tls_genes <- tls_genes[!is.na(tls_genes)]
tls_gene_list <- list(tls_genes)
Idents(smi) <- "tissue"
smi_T  <- subset(smi, idents = "T")
smi <- AddModuleScore(
    object = smi,
    features = list(TLS = tls_genes),
    name = "TLS_Score"
)

## Tumour region
df_tls <- smi@meta.data %>%
    dplyr::select(TLS = TLS_Score1, tissue, group)
df_tls <- df_tls %>% dplyr::filter(tissue %in% c("T"))
df_tls$group <- factor(df_tls$group, levels = c("S", "P", "MP"))
kw_tls <- df_tls %>% kruskal_test(TLS ~ group)
print(kw_tls)
wilcox_tls <- df_tls %>%
    pairwise_wilcox_test(
        TLS ~ group,
        p.adjust.method = "BH"
    ) %>%
  dplyr::rename(
        xmin = group1,
        xmax = group2
    ) %>%
  dplyr::mutate(
      group1 = factor(xmin, levels = levels(df_tls$group)),
        group2 = factor(xmax, levels = levels(df_tls$group)),
        y.position = max(df_tls$TLS) * c(1.05, 1.10, 1.15))
print(wilcox_tls)
nature_colors <- c(
  "S"  = "#4DBBD5FF",  # cyan-blue
  "MP" = "#E64B35FF",  # red
  "P"  = "#00A087FF"   # green
)
comparisons <- list(c("S", "P"), c("S", "MP"), c("P", "MP"))
p_tls <- ggboxplot(df_tls, x = "group", y = "TLS", color = "group", palette = nature_colors, add = "jitter") +
  stat_compare_means(comparisons = comparisons, label = "p.signif", method = "wilcox", p.adjust.method = "BH") +
  stat_compare_means(label.y = 4) +
  scale_fill_manual(values = nature_colors) +
  labs(
    title = "TLS Score Across Groups",
    x = "",
    y = "TLS Score"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_rect(color = "black", size = 1),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )
ggsave(filename = file.path(PROJECT_ROOT, "figures", "TLS_score_tumor_core.pdf"), width = 5, height = 6, plot = p_tls, bg = "white")
