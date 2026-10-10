# -----------------------------------------------------------------------------
# Visium HD - spatial maps of lymphocyte subtypes and the MARCO+ macrophage overlay
#
# Purpose   : Overlay the T/B lymphocyte subtypes onto the histology of every
#             section, and, in a second block, overlay those same lymphocyte
#             subtypes together with the ma_MARCO macrophage subset, to show
#             which myeloid compartment co-localises with the lymphocyte-rich
#             regions.
#
# Inputs    : $PROJECT_ROOT/visium_hd/raw_data/subtype_reannotation.rData
#
# Outputs   : $PROJECT_ROOT/visium_hd/figures/ma_TB_spatial_plots/*.pdf|png
#             $PROJECT_ROOT/visium_hd/figures/ma_TB+MARCO_spatial_plots/*.pdf|png
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 05_spatial_subtype_maps.R
# -----------------------------------------------------------------------------

library(Seurat)
library(ggplot2)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")
fig_dir      <- file.path(HD_ROOT, "figures")

load(file.path(HD_ROOT, "raw_data", "subtype_reannotation.rData"))
# ensure the metadata columns are character to avoid factor-matching issues
object$cell_subtype <- as.character(object$cell_subtype)
object$orig.ident   <- as.character(object$orig.ident)

# =============================================================================
# T/B lymphocyte subtypes (no macrophage overlay)
# =============================================================================
out_dir <- file.path(fig_dir, "ma_TB_spatial_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
target_cells <- c("T_eff","Tfh","T_naive","B_naive","B_active","B_mem","plasma")
TB_colors <- c(
  "T_eff"   = "#E64A3C",   # bright red
  "Tfh"     = "#F6A65B",   # bright orange
  "T_naive" = "#FF6F91",   # bright pink
  "B_naive" = "#A6CEE3",   # pale blue
  "B_active"= "#1F78B4",   # cool blue
  "B_mem"   = "#5E89C7",   # grey-blue
  "plasma"  = "#C2A5CF"    # soft purple-grey
)
# light grey for non-target cells
TB_colors_with_other <- c(
  "Other" = "#EDEDEB",  # neutral background grey
  TB_colors
)
object$TB_group <- ifelse(object$cell_subtype %in% target_cells,
                          object$cell_subtype,
                          "Other")
levels_use <- c("Other", target_cells)
object$TB_group <- factor(object$TB_group, levels = levels_use)
for (img_name in names(object@images)) {
  message("Processing: ", img_name)

  p <- SpatialDimPlot(
    object,
    images = img_name,
    crop = TRUE,
    group.by = "TB_group",
    cols = TB_colors_with_other,
    pt.size.factor = 2,
    alpha = 0.9,
    image.alpha = 0.5
  ) +
    theme_void() +
    ggtitle(paste0(img_name, " - T/B cell subtypes")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.key = element_rect(fill = NA, color = NA)
    )
  fname_base <- file.path(out_dir, paste0(img_name, "_TB_subtypes_spatial"))
  ggsave(paste0(fname_base, ".png"), plot = p, width = 6, height = 6, dpi = 600)
  ggsave(paste0(fname_base, ".pdf"), plot = p, width = 6, height = 6)
}

# =============================================================================
# T/B lymphocyte subtypes + ma_MARCO
# =============================================================================
out_dir <- file.path(fig_dir, "ma_TB+MARCO_spatial_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
target_cells <- c("T_eff","Tfh","T_naive","B_naive","B_active","B_mem","plasma","ma_MARCO")
TB_colors <- c(
  "T_eff"   = "#E64A3C",
  "Tfh"     = "#F6A65B",
  "T_naive" = "#FF6F91",
  "B_naive" = "#A6CEE3",
  "B_active"= "#1F78B4",
  "B_mem"   = "#5E89C7",
  "plasma"  = "#C2A5CF",
  "ma_MARCO" = "#F4D03F"
)
TB_colors_with_other <- c(
  "Other" = "#EDEDEB",
  TB_colors
)
object$TB_group <- ifelse(object$cell_subtype %in% target_cells,
                          object$cell_subtype,
                          "Other")
levels_use <- c("Other", target_cells)
object$TB_group <- factor(object$TB_group, levels = levels_use)
for (img_name in names(object@images)) {
  message("Processing: ", img_name)

  p <- SpatialDimPlot(
    object,
    images = img_name,
    crop = TRUE,
    group.by = "TB_group",
    cols = TB_colors_with_other,
    pt.size.factor = 2,
    alpha = 0.9,
    image.alpha = 0.5
  ) +
    theme_void() +
    ggtitle(paste0(img_name, " - T/B cell subtypes")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.key = element_rect(fill = NA, color = NA)
    )
  fname_base <- file.path(out_dir, paste0(img_name, "_TB+MARCO_subtypes_spatial"))
  ggsave(paste0(fname_base, ".png"), plot = p, width = 6, height = 6, dpi = 600)
  ggsave(paste0(fname_base, ".pdf"), plot = p, width = 6, height = 6)
}
