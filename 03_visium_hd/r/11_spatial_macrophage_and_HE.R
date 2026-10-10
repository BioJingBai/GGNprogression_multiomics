# -----------------------------------------------------------------------------
# Visium HD - spatial localisation of the MARCO+ macrophage subset and raw H&E export
#
# Purpose   : Plot the ma_MARCO macrophage subset alone on top of the histology
#             of every section, and export the raw H&E images for the panel
#             layout.
#
# Inputs    : $PROJECT_ROOT/visium_hd/raw_data/subtype_reannotation.rData
#
# Outputs   : $PROJECT_ROOT/visium_hd/figures/ma_MARCO_spatial_plots/*.pdf|png
#             $PROJECT_ROOT/visium_hd/figures/HE_spatial_plots/*_HE_raw.png
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 11_spatial_macrophage_and_HE.R
# -----------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(png)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")
fig_dir      <- file.path(HD_ROOT, "figures")

load(file.path(HD_ROOT, "raw_data", "subtype_reannotation.rData"))
# ensure the metadata columns are character to avoid factor-matching issues
object$cell_subtype <- as.character(object$cell_subtype)
object$orig.ident   <- as.character(object$orig.ident)

# ---- ma_MARCO ---------------------------------------------------------------
out_dir <- file.path(fig_dir, "ma_MARCO_spatial_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
target_cells <- c("ma_MARCO")
TB_colors <- c("ma_MARCO" = "#F4D03F")
TB_colors_with_other <- c("Other" = "#EDEDEB", TB_colors)
object$TB_group <- ifelse(object$cell_subtype %in% target_cells, object$cell_subtype, "Other")
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
    image.alpha = 0.1
  ) +
    theme_void() +
    ggtitle(paste0(img_name, " - MARCO cell subtypes")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.title = element_blank(),
      legend.key = element_rect(fill = NA, color = NA)
    )
  fname_base <- file.path(out_dir, paste0(img_name, "_MARCO_subtypes_spatial"))
  ggsave(paste0(fname_base, ".png"), plot = p, width = 6, height = 6, dpi = 600)
  ggsave(paste0(fname_base, ".pdf"), plot = p, width = 6, height = 6)
}

# ---- raw H&E export ---------------------------------------------------------
out_dir <- file.path(fig_dir, "HE_spatial_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (img_name in names(object@images)) {
  # raw RGB image
  img <- object@images[[img_name]]@image
  # output filename
  fname <- file.path(out_dir, paste0(img_name, "_HE_raw.png"))
  # write PNG
  writePNG(img, target = fname)
}
for (img_name in names(object@images)) {
 p <- SpatialDimPlot(
  object,
  images = img_name,
  crop = TRUE,
  group.by = NULL,       # no grouping
  pt.size.factor = 0,    # no points
  image.alpha = 1        # original brightness
) + theme_void() +
    ggtitle(img_name)
  fname_base <- file.path(out_dir, paste0(img_name, "_HE"))
  png(filename = fname, width = 6000, height = 6000, res = 600)  # very high resolution
  print(p)  # SpatialDimPlot object
  dev.off()
  ggsave(paste0(fname_base, ".pdf"), plot = p, width = 6, height = 6)
}
