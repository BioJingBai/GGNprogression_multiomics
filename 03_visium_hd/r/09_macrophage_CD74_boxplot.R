# -----------------------------------------------------------------------------
# Visium HD - CD74 expression across macrophage subsets
#
# Purpose   : Compare CD74 expression between ma_MARCO and the other macrophage
#             subsets (ma_NR4A1, ma_SPP1, ma_FOLR2, ma_EPAS1) with Wilcoxon
#             rank-sum tests (BH-adjusted) and visualise the distribution per
#             subset.
#
# Inputs    : $PROJECT_ROOT/visium_hd/raw_data/subtype_reannotation.rData
#
# Outputs   : $PROJECT_ROOT/visium_hd/figures/CD74 in ma_marco.pdf
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 09_macrophage_CD74_boxplot.R
# -----------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(FSA)      # dunnTest
library(readr)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")
dir.create(file.path(HD_ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)

load(file.path(HD_ROOT, "raw_data", "subtype_reannotation.rData"))
All_Seurat <- object
macro_types <- c("ma_EPAS1", "ma_FOLR2", "ma_SPP1", "ma_MARCO", "ma_NR4A1")
df <- FetchData(All_Seurat, vars = c("CD74", "cell_subtype")) |>
  filter(cell_subtype %in% macro_types)
dunn_res <- dunnTest(`CD74` ~ cell_subtype, data = df, method = "bonferroni")
dunn_res
# comparison pairs to annotate
comparisons <- list(
  c("ma_MARCO", "ma_NR4A1"),
  c("ma_MARCO", "ma_SPP1"),
  c("ma_MARCO", "ma_FOLR2"),
  c("ma_MARCO", "ma_EPAS1")
)
pdf(file.path(HD_ROOT, "figures", "CD74 in ma_marco.pdf"), width = 6, height = 6)
ggplot(df, aes(x = cell_subtype, y = `CD74`, fill = cell_subtype)) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  theme_classic() +
  labs(x = "",
       y = "CD74 expression",
       title = "CD74 distribution across macrophage subclusters") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
