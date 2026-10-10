# -----------------------------------------------------------------------------
# Distribution of SPP1+ macrophages across regions and groups
#
# Purpose   : Count SPP1+ macrophages (final_Level4_celltype == "ma_SPP1") per
#             Group7 (group x region) and per sample, plot the bar charts and
#             run the pairwise Wilcoxon/t-test comparison between groups.
#
# Inputs    : <PROJECT_ROOT>/processed/smi_annotated.qs
#
# Outputs   : <PROJECT_ROOT>/figures/SPP1mac_compare.pdf
#             <PROJECT_ROOT>/figures/SPP1mac_compare_per_sample.pdf
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 03_figure6c_spp1_macrophage_distribution.R
# -----------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Dependencies
## ---------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(qs)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
dir.create(file.path(PROJECT_ROOT, "figures"), recursive = TRUE, showWarnings = FALSE)

## Load the annotated object produced by 00_setup_and_annotation.R
smi <- qread(file.path(PROJECT_ROOT, "processed", "smi_annotated.qs"))


## ---------------------------------------------------------------------------
## SPP1+ macrophages in each region
## ---------------------------------------------------------------------------
mac_cells <- rownames(smi@meta.data)[smi$final_Level2_celltype == "Macrophages"]
## Summary over all samples
df <- smi@meta.data %>%
  dplyr::filter(smi$final_Level4_celltype == "ma_SPP1") %>%
  group_by(Group7) %>%
  dplyr::summarise(n = n())
pdf(file.path(PROJECT_ROOT, "figures", "SPP1mac_compare.pdf"), height = 6, width = 8)
ggplot(df, aes(x = Group7, y = n, fill = Group7)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  ylab("Number of SPP1+ Macrophages") +
  xlab("Group") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ggtitle("Distribution of SPP1+ Macrophages across groups")
dev.off()

## Compare between groups per sample (not significant)
df <- smi@meta.data %>%
  filter(smi$final_Level4_celltype == "ma_SPP1") %>%
  group_by(Group7, sample) %>%   # split by group and sample
  dplyr::summarise(n = n()) %>%
  ungroup()
pdf(file.path(PROJECT_ROOT, "figures", "SPP1mac_compare_per_sample.pdf"), height = 6, width = 8)
# Bar chart: each bar is the within-group mean
ggplot(df, aes(x = Group7, y = n, fill = Group7)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge()) +  # bars are the group means
  geom_jitter(aes(color = Group7), width = 0.2, size = 2) +  # sample points
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  ylab("Number of SPP1+ Macrophages per sample") +
  xlab("Group") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ggtitle("Distribution of SPP1+ Macrophages per sample")
dev.off()

# Between-group t-test or non-parametric test
# e.g. all pairwise comparisons:
pairwise_results <- df %>%
  rstatix::pairwise_t_test(
    n ~ Group7,         # compare the counts between groups
    p.adjust.method = "BH"  # multiple-testing correction
  )
# Inspect the result
pairwise_results
