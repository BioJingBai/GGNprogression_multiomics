# -----------------------------------------------------------------------------
# Visium HD - tissue preference (Ro/e) of major cell types across groups
#
# Purpose   : Compute the observed-over-expected ratio (Ro/e) of every major
#             cell type across the three progression groups and assemble the
#             composite figure: a stacked top bar (group composition), a Ro/e
#             heatmap in the middle and a stacked right bar (overall cell-type
#             composition).
#
# Inputs    : $PROJECT_ROOT/visium_hd/raw_data/subtype_reannotation.rData
#
# Outputs   : $PROJECT_ROOT/visium_hd/figures/celltypr barplot.pdf
#             $PROJECT_ROOT/visium_hd/figures/Ro_e_visium_hd_major_celltypes.csv
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 04_roie_tissue_preference_heatmap.R
# -----------------------------------------------------------------------------

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")
dir.create(file.path(HD_ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)

load(file.path(HD_ROOT, "raw_data", "subtype_reannotation.rData"))
meta <- object@meta.data
table(meta$cell_type_renew_Retain_malignant, meta$Group)
meta <- dplyr::filter(meta, meta$Group != "S0432519_9")
prop_df <- meta %>%
  group_by(Group, cell_type_renew_Retain_malignant) %>%
  summarise(n = n()) %>%
  group_by(Group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
prop_df$group <- ifelse(prop_df$Group == "S0439254_9", "MP",
                        ifelse(prop_df$Group == "S0441838_5", "P", "S"))
prop_mat <- prop_df %>%
  dplyr::select(cell_type_renew_Retain_malignant, group, prop) %>%
  pivot_wider(names_from = group, values_from = prop, values_fill = 0)
prop_mat <- as.data.frame(prop_mat)
rownames(prop_mat) <- prop_mat$cell_type_renew_Retain_malignant
prop_mat$cell_type_renew_Retain_malignant <- NULL
prop_mat <- as.matrix(prop_mat)[, c("MP", "P", "S")]

# ---- Ro/e -------------------------------------------------------------------
# The same cell-type x group table, now in observed counts.
count_mat <- prop_df %>%
  dplyr::select(cell_type_renew_Retain_malignant, group, n) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0)
count_mat <- as.data.frame(count_mat)
rownames(count_mat) <- count_mat$cell_type_renew_Retain_malignant
count_mat$cell_type_renew_Retain_malignant <- NULL
count_mat <- as.matrix(count_mat)[, c("MP", "P", "S")]

# Expected counts under the null hypothesis that cell type and group are
# independent, taken from the chi-square test of the contingency table
# (Methods: "the expected counts under the null hypothesis were calculated
# using a chi-square test").
chisq_res <- suppressWarnings(chisq.test(count_mat))
exp_mat   <- chisq_res$expected

# Ro/e = observed / expected: >1 enriched in that group, <1 depleted.
Ro_e <- count_mat / exp_mat

write.csv(
  data.frame(cell_type = rownames(Ro_e), round(Ro_e, 4), check.names = FALSE),
  file.path(HD_ROOT, "figures", "Ro_e_visium_hd_major_celltypes.csv"),
  row.names = FALSE
)

celltype_colors <- c(
  "B_cell"="#66C2A5",
  "CAFs"="#FC8D62",
  "cDC"="#8DA0CB",
  "Endothelial"="#E78AC3",
  "Epithelial"="#A6D854",
  "Macrophage"="#FFD92F",
  "Malignant"="#E5C494",
  "Pericyte"="#B3B3B3",
  "Plasma"="#1F78B4",
  "SMCs"="#FB9A99",
  "T_cell"="#CAB2D6"
)

group_colors <- c(
  "MP"="#C27AA8",
  "P"="#3970A4",
  "S"="#E0D395"
)

# top stacked bar
df_top <- prop_mat %>%
  as.data.frame() %>%
  mutate(celltype = rownames(.)) %>%
  pivot_longer(-celltype, names_to = "group", values_to = "prop")
p_top <- ggplot(df_top, aes(x = group, y = prop, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  ylab("Percentage") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# middle heatmap
row_names <- rownames(Ro_e)
row_ha <- rowAnnotation(
  type = anno_simple(row_names, col = celltype_colors, gp = gpar(col = NA)),
  show_annotation_name = FALSE,
  width = unit(6, "mm")
)
# Diverging scale centred on 1 (white = no enrichment), symmetric on the log2
# axis so that a 2-fold enrichment and a 2-fold depletion match in colour.
ro_lim  <- max(abs(log2(Ro_e)))
col_fun <- colorRamp2(2^c(-ro_lim, 0, ro_lim), c("#424492", "white", "#E28D44"))
ro_at   <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
ro_at   <- ro_at[ro_at >= min(Ro_e) & ro_at <= max(Ro_e)]
if (length(ro_at) == 0) ro_at <- 1
ht <- Heatmap(
  Ro_e,
  name = "R_o/e",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  left_annotation = row_ha,
  heatmap_legend_param = list(
    title = expression(R[o]/e),
    at = ro_at,
    labels = as.character(ro_at),
    legend_direction = "horizontal"
  )
)

# right stacked bar
comp_right <- rowSums(prop_mat)
comp_right <- comp_right / sum(comp_right)
df_right <- data.frame(
  celltype = names(comp_right),
  value = comp_right
)
p_right <- ggplot(df_right, aes(x = 1, y = value, fill = celltype)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = celltype_colors) +
  theme_void() +
  theme(legend.position = "none")

# bottom group colour strip
df_bottom <- data.frame(
  group = names(group_colors),
  color = group_colors
)
p_bottom <- ggplot(df_bottom, aes(x = group, y = 1, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = group_colors) +
  theme_void() +
  theme(legend.position = "none")

# assemble
ht_list <- draw(ht, merge_legend = TRUE, newpage = FALSE)
ht_grob <- grid.grabExpr(draw(ht, merge_legend = TRUE, newpage = FALSE))
ht_gg <- as.ggplot(ht_grob)

p_top <- p_top + plot_layout(guides = "collect")
p_bottom <- p_bottom + plot_layout(guides = "collect")
final_plot <- (
    p_top /
    ht_gg  /
    p_bottom +
    plot_layout(heights = c(1,4,0.4), widths = c(4,1))
)

final_plot
ggsave(file.path(HD_ROOT, "figures", "celltypr barplot.pdf"), final_plot, width = 5, height = 6)
