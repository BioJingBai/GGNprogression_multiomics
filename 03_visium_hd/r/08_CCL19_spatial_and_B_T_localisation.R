# -----------------------------------------------------------------------------
# Visium HD - CCL19 spatial expression and B_naive / T_eff co-localisation
#
# Purpose   : (1) Show CCL19 expression on the Stable section (S0442225_10);
#             (2) highlight B_naive and T_eff cells across all sections, both
#             as a 2x2 combined panel and as one spatial plot per section.
#
# Inputs    : $PROJECT_ROOT/visium_hd/raw_data/subtype_reannotation.rData
#
# Outputs   : CCL19_S0442225_10_Spatial.pdf
#             Bnaive_Teff_spatial_2x2.pdf
#             spatial_<section>.pdf (one per section)
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 08_CCL19_spatial_and_B_T_localisation.R
# -----------------------------------------------------------------------------

library(Seurat)
library(SeuratObject)
library(tidyr)
library(dplyr)
library(Matrix)
library(scales)
library(ggplot2)
library(future)
library(qs)
library(harmony)
library(openxlsx)
library(patchwork)
library(pheatmap)
library(ComplexHeatmap)
library(purrr)

HD_RAW <- file.path(Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data"),
                    "visium_hd", "raw_data")

# ==================== 1. CCL19 expression on the Stable section ====================
load(file.path(HD_RAW, "subtype_reannotation.rData"))
marker_genes <- c("CCL19")
p_spatial <- SpatialFeaturePlot(
  object,
  features = marker_genes[1],
  images = "S0442225_10",
  image.alpha = 0.2,  # background opacity
  pt.size.factor = 2,  # point size
  alpha = c(0.2, 0.6)   # low / high expression point opacity
) +
  scale_fill_gradientn(
    colors = c("white", "#FFFFFF", "#3F3FFF", "#0000FF"),  # same four-colour gradient as the UMAP
    values = scales::rescale(c(0, 0.2, 0.6, 1)),  # colour break points
    limits = c(0, 1),
    name = "Expression\nLevel",
    guide = guide_colorbar(barwidth = 0.8, barheight = 4)  # legend size
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)  # title style
  ) +
  labs(title = "CCL19 Spatial Expression")

# show the image
print(p_spatial)

# save PDF
ggsave("CCL19_S0442225_10_Spatial.pdf", plot = p_spatial, width = 12, height = 5, dpi = 300)

# ==================== 2. B_naive and T_eff distribution ====================
load(file.path(HD_RAW, "subtype_reannotation.rData"))
str(object, 3)
# required packages
library(Seurat)
library(patchwork)   # for laying out multiple ggplots

# 1) add a highlight column in meta.data (B_naive / T_eff / other)
object$highlight <- ifelse(object$cell_subtype == "B_naive", "B_naive",
                           ifelse(object$cell_subtype == "T_eff", "T_eff", "other"))
object$highlight <- factor(object$highlight, levels = c("B_naive", "T_eff", "other"))

# 2) colours (B_naive red, T_eff blue, other grey)
cols <- c("B_naive" = "red", "T_eff" = "blue", "other" = "lightgrey")

# 3) take every sample (section) and make one SpatialDimPlot each
samples <- unique(object$orig.ident)   # alternatively object$Group or the image names, depending on the fields
plots <- lapply(samples, function(samp) {
  subobj <- subset(object, subset = orig.ident == samp)
  # SpatialDimPlot overlays the coordinates on the tissue image
  p <- SpatialDimPlot(subobj,
                      group.by = "highlight",
                      image.alpha = 0.2,  # background opacity
                      cols = cols,
                      # label = FALSE,       # uncomment to hide the class labels
                      repel = TRUE) +       # nicer legend text placement (if supported)
       ggtitle(as.character(samp)) +
       theme(plot.title = element_text(hjust = 0.5))
  return(p)
})

# 4) combine the plots into 2x2 (wrap_plots adapts to the number of samples)
combined <- wrap_plots(plots, ncol = 2) + plot_layout(guides = "collect") &
            theme(legend.position = "bottom")

# 5) save to file
ggsave("Bnaive_Teff_spatial_2x2.pdf", combined, width = 12, height = 10, dpi = 300)

## if 2x2 is not needed, the per-sample layout below can be used instead

# 3) sample list (split by orig.ident)
samples <- unique(object$orig.ident)

# 4) one spatial plot per sample, collected in a list
plots <- lapply(samples, function(samp) {
  subobj <- subset(object, subset = orig.ident == samp)
  p <- SpatialDimPlot(subobj,
                      group.by = "highlight",
                      image.alpha = 0.2,  # background opacity
                      cols = cols,
                      # if points are too dense or too small, adjust pt.size.factor (e.g. 1.2, 0.8)
                      pt.size.factor = 2) +
       ggtitle(as.character(samp)) +
       theme(plot.title = element_text(hjust = 0.5),
             legend.position = "bottom")
  return(p)
})

# 5) print the four plots in sequence
for (p in plots) print(p)

# 6) save each to a file, named by sample
for (i in seq_along(plots)) {
  fname <- paste0("spatial_", gsub("[^A-Za-z0-9_\\-]", "_", samples[i]), ".pdf")
  ggsave(filename = fname, plot = plots[[i]], width = 6, height = 6, dpi = 300)
}
