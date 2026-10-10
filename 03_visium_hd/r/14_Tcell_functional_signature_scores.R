# -----------------------------------------------------------------------------
# Visium HD - functional signature scores across T cell subtypes
#
# Purpose   : Score every T cell for a cytotoxicity signature and an exhaustion
#             signature (AddModuleScore), then compare the two scores between
#             the three progression sections within each T cell subtype
#             (T_eff, Tfh, T_ex). A naive/memory marker dotplot per
#             section-subtype is produced alongside.
#
# Inputs    : $PROJECT_ROOT/visium_hd/raw_data/subtype_reannotation.rData
#
# Outputs   : $PROJECT_ROOT/visium_hd/figures/T_cell_naive_memory_markers.pdf
#             $PROJECT_ROOT/visium_hd/figures/T_eff_cytotoxicity_score.pdf
#             $PROJECT_ROOT/visium_hd/figures/T_eff_exhaustion_score.pdf
#             $PROJECT_ROOT/visium_hd/figures/Tfh_cytotoxicity_score.pdf
#             $PROJECT_ROOT/visium_hd/figures/T_ex_exhaustion_score.pdf
#
# Usage     : PROJECT_ROOT=/path/to/GGN_progression_data Rscript 14_Tcell_functional_signature_scores.R
# -----------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
HD_ROOT      <- file.path(PROJECT_ROOT, "visium_hd")
dir.create(file.path(HD_ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)

load(file.path(HD_ROOT, "raw_data", "subtype_reannotation.rData"))

samples_to_use <- c("S0439254_9","S0441838_5", "S0442225_10")
T_cells <- subset(object, subset = cell_type_renew == "T_cell" & orig.ident %in% samples_to_use)
T_cells$sample_Tsub <- paste(T_cells$orig.ident, T_cells$cell_subtype, sep = "_")
T_cell_naive_mem_markers <- c("SELL","TCF7","CCR7","IL7R","ANXA1","GPR183","GZMK","GZMB","PRF1","NKG7","GZMH","GNLY","XCL1","XCL2","LAG3","HAVCR2","TNFRSF9","TIGIT","PDCD1","ENTPD1","TOX")
p <- DotPlot(T_cells, features = T_cell_naive_mem_markers, group.by = "sample_Tsub", scale.by = "size", scale = 10, col.max = 1.5) +
        RotatedAxis() +
        scale_color_gradientn(values = seq(0,1,0.1), colours = c("#4575b4","#abd9e9","#e0f3f8","#ffffbf","#fdae61","#d73027","#800026"))
ggsave(file.path(HD_ROOT, "figures", "T_cell_naive_memory_markers.pdf"), p, width = 12, height = 5)


T_cytotoxicity_signature <- list(c("GZMK", "GPR183", "GNLY", "GZMH", "CX3CR1", "NKG7", "TBX21", "PRF1", "CD300A", "EOMES"))
T_exhaustion_signature <- list(c("GZMB", "HAVCR2", "ENTPD1", "VCAM1", "LAG3", "TNFRSF9", "CTLA4", "TIGIT", "TNFRSF18", "PDCD1", "IFNG", "CXCR6", "GNLY", "DUSP4"))
T_cells <- AddModuleScore(T_cells, features = T_cytotoxicity_signature, name = "T_cytotoxicity_signature", nbin = 10, ctrl = 20)
T_cells <- AddModuleScore(T_cells, features = T_exhaustion_signature, name = "T_exhaustion_signature", nbin = 10, ctrl = 20)
T_cell_exp <- as.data.frame(T_cells@meta.data)[, c(1,14,18,19)]
# define functions for half violin plot #
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                        draw_group = function(self, data, ..., draw_quantiles = NULL) {
                                data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                                grp <- data[1, "group"]
                                newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                                newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                                newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                                if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                            1))
                                quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                aesthetics$alpha <- rep(1, nrow(quantiles))
                                both <- cbind(quantiles, aesthetics)
                                quantile_grob <- GeomPath$draw_panel(both, ...)
                                ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                                }
                                else {
                                ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                                }
                            }
                        )
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                            draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                            show.legend = NA, inherit.aes = TRUE) {
                                                                    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
                                                                    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                                                                    params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
                                                                    }

# plot #
T_cell_exp_eff <- dplyr::filter(T_cell_exp, T_cell_exp$cell_subtype == "T_eff")
comparisons <- list(c("S0442225_10","S0441838_5"), c("S0442225_10","S0439254_9"), c("S0441838_5","S0439254_9"))
p <- ggplot(T_cell_exp_eff, aes(x = orig.ident, y = T_cytotoxicity_signature1, fill = orig.ident)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.15, outliers = F, fill = "white") +
    scale_fill_manual(values = c("#a6cee3","#fdbf6f","#4a8594")) +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  stat_compare_means(comparisons = comparisons, label = "p.signif", method = "wilcox", p.adjust.method = "BH") +
      stat_compare_means(label.y = 0.28) +
    ggtitle('T_cytotoxicity_signature') +
    guides(fill = guide_legend(title = NULL)) +
    xlab("group") +
    ylab("Signature score")
p <- ggsave(file.path(HD_ROOT, "figures", "T_eff_cytotoxicity_score.pdf"), p, width = 4, height = 6)
T_cell_exp_eff %>%
  group_by(orig.ident) %>%
  summarise(
    n = n(),
    mean_score = mean(T_cytotoxicity_signature1, na.rm = TRUE),
    median_score = median(T_cytotoxicity_signature1, na.rm = TRUE),
    sd_score = sd(T_cytotoxicity_signature1, na.rm = TRUE),
    IQR_score = IQR(T_cytotoxicity_signature1, na.rm = TRUE),
    min_score = min(T_cytotoxicity_signature1, na.rm = TRUE),
    max_score = max(T_cytotoxicity_signature1, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_score))
pairwise.t.test(
  x = T_cell_exp_eff$T_cytotoxicity_signature1,
  g = T_cell_exp_eff$orig.ident,
  p.adjust.method = "BH"
)
# plot #
T_cell_exp_eff <- dplyr::filter(T_cell_exp, T_cell_exp$cell_subtype == "T_eff")
comparisons <- list(c("S0442225_10","S0441838_5"), c("S0442225_10","S0439254_9"), c("S0441838_5","S0439254_9"))
p <- ggplot(T_cell_exp_eff, aes(x = orig.ident, y = T_exhaustion_signature1, fill = orig.ident)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.15, outliers = F, fill = "white") +
    scale_fill_manual(values = c("#a6cee3","#fdbf6f","#4a8594")) +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  stat_compare_means(comparisons = comparisons, label = "p.signif", method = "wilcox", p.adjust.method = "BH") +
      stat_compare_means(label.y = 0.28) +
    ggtitle('T_exhaustion_signature') +
    guides(fill = guide_legend(title = NULL)) +
    xlab("group") +
    ylab("Signature score")
p <- ggsave(file.path(HD_ROOT, "figures", "T_eff_exhaustion_score.pdf"), p, width = 4, height = 6)
T_cell_exp_eff %>%
  group_by(orig.ident) %>%
  summarise(
    n = n(),
    mean_score = mean(T_exhaustion_signature1, na.rm = TRUE),
    median_score = median(T_exhaustion_signature1, na.rm = TRUE),
    sd_score = sd(T_exhaustion_signature1, na.rm = TRUE),
    IQR_score = IQR(T_exhaustion_signature1, na.rm = TRUE),
    min_score = min(T_exhaustion_signature1, na.rm = TRUE),
    max_score = max(T_exhaustion_signature1, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_score))
pairwise.t.test(
  x = T_cell_exp_eff$T_exhaustion_signature1,
  g = T_cell_exp_eff$orig.ident,
  p.adjust.method = "BH"
)

# plot #
T_cell_exp_tfh <- dplyr::filter(T_cell_exp, T_cell_exp$cell_subtype == "Tfh")
comparisons <- list(c("S0442225_10","S0441838_5"), c("S0442225_10","S0439254_9"), c("S0441838_5","S0439254_9"))
p <- ggplot(T_cell_exp_tfh, aes(x = orig.ident, y = T_exhaustion_signature1, fill = orig.ident)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.15, outliers = F, fill = "white") +
    scale_fill_manual(values = c("#a6cee3","#fdbf6f","#4a8594")) +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  stat_compare_means(comparisons = comparisons, label = "p.signif", method = "wilcox", p.adjust.method = "BH") +
      stat_compare_means(label.y = 0.28) +
    ggtitle('T_exhaustion_signature') +
    guides(fill = guide_legend(title = NULL)) +
    xlab("group") +
    ylab("Signature score")
p <- ggsave(file.path(HD_ROOT, "figures", "Tfh_cytotoxicity_score.pdf"), p, width = 4, height = 6)
T_cell_exp_tfh %>%
  group_by(orig.ident) %>%
  summarise(
    n = n(),
    mean_score = mean(T_cytotoxicity_signature1, na.rm = TRUE),
    median_score = median(T_cytotoxicity_signature1, na.rm = TRUE),
    sd_score = sd(T_cytotoxicity_signature1, na.rm = TRUE),
    IQR_score = IQR(T_cytotoxicity_signature1, na.rm = TRUE),
    min_score = min(T_cytotoxicity_signature1, na.rm = TRUE),
    max_score = max(T_cytotoxicity_signature1, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_score))

pairwise.t.test(
  x = T_cell_exp_tfh$T_cytotoxicity_signature1,
  g = T_cell_exp_tfh$orig.ident,
  p.adjust.method = "BH"
)

# plot #
T_cell_exp_ex <- dplyr::filter(T_cell_exp, T_cell_exp$cell_subtype == "T_ex")
comparisons <- list(c("S0442225_10","S0441838_5"), c("S0442225_10","S0439254_9"), c("S0441838_5","S0439254_9"))
p <- ggplot(T_cell_exp_ex, aes(x = orig.ident, y = T_exhaustion_signature1, fill = orig.ident)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.15, outliers = F, fill = "white") +
    scale_fill_manual(values = c("#a6cee3","#fdbf6f","#4a8594")) +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  stat_compare_means(comparisons = comparisons, label = "p.signif", method = "wilcox", p.adjust.method = "BH") +
      stat_compare_means(label.y = 0.18) +
    ggtitle('T_exhaustion_signature') +
    guides(fill = guide_legend(title = NULL)) +
    xlab("group") +
    ylab("Signature score")
p <- ggsave(file.path(HD_ROOT, "figures", "T_ex_exhaustion_score.pdf"), p, width = 4, height = 6)
T_cell_exp_ex %>%
  group_by(orig.ident) %>%
  summarise(
    n = n(),
    mean_score = mean(T_exhaustion_signature1, na.rm = TRUE),
    median_score = median(T_exhaustion_signature1, na.rm = TRUE),
    sd_score = sd(T_exhaustion_signature1, na.rm = TRUE),
    IQR_score = IQR(T_exhaustion_signature1, na.rm = TRUE),
    min_score = min(T_exhaustion_signature1, na.rm = TRUE),
    max_score = max(T_exhaustion_signature1, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_score))


pairwise.t.test(
  x = T_cell_exp_ex$T_cytotoxicity_signature1,
  g = T_cell_exp_ex$orig.ident,
  p.adjust.method = "BH"
)
