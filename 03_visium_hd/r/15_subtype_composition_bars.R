# -----------------------------------------------------------------------------
# Visium HD - subtype composition stacked bars per section
#
# Purpose   : Turn the per-section subtype counts into proportion stacked bars
#             (alluvial-style, one bar per section) for the B, macrophage, T and
#             epithelial compartments.
#
# Inputs    : $PROJECT_ROOT/visium_hd/subtype_counts/{B_cell_subtypes,Macrophage_subtypes,T_cell_subtypes,epithelial_subtypes}/*_subtypes_counts_by_sample.csv
#
# Outputs   : $PROJECT_ROOT/visium_hd/figures/subtype_composition/*_subtypes_stacked_bar_plot.pdf|png
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 15_subtype_composition_bars.R
# -----------------------------------------------------------------------------

rm(list = ls())
gc()

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(tidyverse)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")
out_dir      <- file.path(HD_ROOT, "figures", "subtype_composition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# data directory
base_dir <- file.path(HD_ROOT, "subtype_counts")

# colour scheme
# palette 5 - custom colours (9 colours)
mycolor <- c('#E59572', '#E7D1C9', '#FF8EAD', '#FFEEAD', '#D9534F',
             '#FFAD60', '#57D1C9', '#F4A460', '#98D8C8')

# ==================== 1. B cell subtypes ====================
message("Processing: B cell subtypes")

# read data
b_data <- read.csv(file.path(base_dir, "B_cell_subtypes", "b_subtypes_counts_by_sample.csv"),
                   header = TRUE, sep = ",")

# data wrangling
b_data_ratio <- b_data %>%
  mutate(total = rowSums(across(-orig.ident))) %>%
  mutate(across(-c(orig.ident, total), ~ ./total)) %>%
  select(-total) %>%
  pivot_longer(cols = -orig.ident, names_to = "celltype", values_to = "proportion") %>%
  pivot_wider(names_from = orig.ident, values_from = "proportion") %>%
  pivot_longer(-celltype, names_to = "group", values_to = "value") %>%
  mutate(group = factor(group, levels = c("S0432519_9", "S0439254_9", "S0441838_5", "S0442225_10")))

# plot
p_b <- ggplot(b_data_ratio, aes(x = group, y = value, fill = celltype,
                                stratum = celltype, alluvium = celltype)) +
  geom_col(width = 0.5, color = 'black') +
  geom_flow(width = 0.5, alpha = 0.4, knot.pos = 0.5) +
  theme_classic() +
  labs(x = '', y = 'proportion', title = "") +
  scale_fill_manual(values = mycolor[1:length(unique(b_data_ratio$celltype))]) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# save
ggsave(file.path(out_dir, "b_subtypes_stacked_bar_plot.pdf"), p_b, width = 4, height = 6)
ggsave(file.path(out_dir, "b_subtypes_stacked_bar_plot.png"), p_b, width = 4, height = 6, dpi = 300)
message("B cell subtypes done")

# ==================== 2. Macrophage subtypes ====================
message("Processing: Macrophage subtypes")

# read data
macro_data <- read.csv(file.path(base_dir, "Macrophage_subtypes", "Macrophage_subtypes_counts_by_sample.csv"),
                       header = TRUE, sep = ",")

# data wrangling
macro_data_ratio <- macro_data %>%
  mutate(total = rowSums(across(-orig.ident))) %>%
  mutate(across(-c(orig.ident, total), ~ ./total)) %>%
  select(-total) %>%
  pivot_longer(cols = -orig.ident, names_to = "celltype", values_to = "proportion") %>%
  pivot_wider(names_from = orig.ident, values_from = "proportion") %>%
  pivot_longer(-celltype, names_to = "group", values_to = "value") %>%
  mutate(group = factor(group, levels = c("S0432519_9", "S0439254_9", "S0441838_5", "S0442225_10")))

# plot
p_macro <- ggplot(macro_data_ratio, aes(x = group, y = value, fill = celltype,
                                        stratum = celltype, alluvium = celltype)) +
  geom_col(width = 0.5, color = 'black') +
  geom_flow(width = 0.5, alpha = 0.4, knot.pos = 0.5) +
  theme_classic() +
  labs(x = '', y = 'proportion', title = "") +
  scale_fill_manual(values = mycolor[1:length(unique(macro_data_ratio$celltype))]) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# save
ggsave(file.path(out_dir, "Macrophage_subtypes_stacked_bar_plot.pdf"), p_macro, width = 4, height = 6)
ggsave(file.path(out_dir, "Macrophage_subtypes_stacked_bar_plot.png"), p_macro, width = 4, height = 6, dpi = 300)
message("Macrophage subtypes done")

# ==================== 3. T cell subtypes ====================
message("Processing: T cell subtypes")

# read data
t_data <- read.csv(file.path(base_dir, "T_cell_subtypes", "T_subtypes_counts_by_sample.csv"),
                   header = TRUE, sep = ",")

# data wrangling
t_data_ratio <- t_data %>%
  mutate(total = rowSums(across(-orig.ident))) %>%
  mutate(across(-c(orig.ident, total), ~ ./total)) %>%
  select(-total) %>%
  pivot_longer(cols = -orig.ident, names_to = "celltype", values_to = "proportion") %>%
  pivot_wider(names_from = orig.ident, values_from = "proportion") %>%
  pivot_longer(-celltype, names_to = "group", values_to = "value") %>%
  mutate(group = factor(group, levels = c("S0432519_9", "S0439254_9", "S0441838_5", "S0442225_10")))

# plot
p_t <- ggplot(t_data_ratio, aes(x = group, y = value, fill = celltype,
                                stratum = celltype, alluvium = celltype)) +
  geom_col(width = 0.5, color = 'black') +
  geom_flow(width = 0.5, alpha = 0.4, knot.pos = 0.5) +
  theme_classic() +
  labs(x = '', y = 'proportion', title = "") +
  scale_fill_manual(values = mycolor[1:length(unique(t_data_ratio$celltype))]) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# save
ggsave(file.path(out_dir, "T_subtypes_stacked_bar_plot.pdf"), p_t, width = 5, height = 6)
ggsave(file.path(out_dir, "T_subtypes_stacked_bar_plot.png"), p_t, width = 5, height = 6, dpi = 300)
message("T cell subtypes done")

# ==================== 4. Epithelial subtypes ====================
message("Processing: epithelial subtypes")

# read data
epi_data <- read.csv(file.path(base_dir, "epithelial_subtypes", "Epi_subtypes_counts_by_sample.csv"),
                     header = TRUE, sep = ",")

# data wrangling
epi_data_ratio <- epi_data %>%
  mutate(total = rowSums(across(-orig.ident))) %>%
  mutate(across(-c(orig.ident, total), ~ ./total)) %>%
  select(-total) %>%
  pivot_longer(cols = -orig.ident, names_to = "celltype", values_to = "proportion") %>%
  pivot_wider(names_from = orig.ident, values_from = "proportion") %>%
  pivot_longer(-celltype, names_to = "group", values_to = "value") %>%
  mutate(group = factor(group, levels = c("S0432519_9", "S0439254_9", "S0441838_5", "S0442225_10")))

# plot
p_epi <- ggplot(epi_data_ratio, aes(x = group, y = value, fill = celltype,
                                    stratum = celltype, alluvium = celltype)) +
  geom_col(width = 0.5, color = 'black') +
  geom_flow(width = 0.5, alpha = 0.4, knot.pos = 0.5) +
  theme_classic() +
  labs(x = '', y = 'proportion', title = "") +
  scale_fill_manual(values = mycolor[1:length(unique(epi_data_ratio$celltype))]) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# save
ggsave(file.path(out_dir, "Epi_subtypes_stacked_bar_plot.pdf"), p_epi, width = 5, height = 6)
ggsave(file.path(out_dir, "Epi_subtypes_stacked_bar_plot.png"), p_epi, width = 5, height = 6, dpi = 300)
message("Epithelial subtypes done")

message("All subtypes processed!")
