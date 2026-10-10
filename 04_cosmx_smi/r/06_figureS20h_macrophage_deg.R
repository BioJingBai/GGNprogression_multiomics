# -----------------------------------------------------------------------------
# Differential gene expression of tumour-core macrophages
#
# Purpose   : Compare the tumour-core (tissue == "T") macrophages between the MP
#             and P groups (FindMarkers, logistic regression) and draw the
#             volcano plot with the top up- and down-regulated genes labelled.
#
# Inputs    : <PROJECT_ROOT>/processed/smi_annotated.qs
#
# Outputs   : <PROJECT_ROOT>/figures/macrophage_DEG_tumor_core_MP_vs_P.pdf
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 06_figureS20h_macrophage_deg.R
# -----------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Dependencies
## ---------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(qs)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
dir.create(file.path(PROJECT_ROOT, "figures"), recursive = TRUE, showWarnings = FALSE)

## Load the annotated object produced by 00_setup_and_annotation.R
smi <- qread(file.path(PROJECT_ROOT, "processed", "smi_annotated.qs"))


## ---------------------------------------------------------------------------
## Differential expression: tumour-core macrophages, MP vs P
## ---------------------------------------------------------------------------
# Level2_celltype lives in meta.data and Sample in meta.data as well
macrophages_T <- subset(smi, subset = Level2_celltype == "Macrophages" & tissue == "T")
DE_genes <- FindMarkers(
  object = macrophages_T,
  ident.1 = "MP",
  ident.2 = "P",
  group.by = "group",   # specify the grouping column
  logfc.threshold = 0.25,  # threshold can be tuned
  min.pct = 0.1,           # expressed in at least 10% of the cells
  test.use = "LR"
)
DE_genes$gene <- rownames(DE_genes)
diff_g1 <- DE_genes
diff_g1$g <- "normal"
diff_g1$g[diff_g1$avg_log2FC > 0 & diff_g1$p_val_adj < 0.05] <- "up"
diff_g1$g[diff_g1$avg_log2FC < 0 & diff_g1$p_val_adj < 0.05] <- "down"
head(diff_g1)
table(diff_g1$g)
diff_g1$color <- paste0(diff_g1$g, diff_g1$targets)
table(diff_g1$color)
# Extract the non-significant and the significant points
data_bg <- diff_g1[diff_g1$g == "normal", ]
data_sig <- diff_g1[diff_g1$g != "normal", ]
data_sig$color <- factor(data_sig$color, levels = c("up", "down"))
# Top genes to label
# Significant differentially expressed genes
top10 <- data_sig %>%
  arrange(desc(avg_log2FC)) %>%   # could also use desc(logFC)
  slice(1:5)
bottom10 <- data_sig %>%
  arrange(avg_log2FC) %>%   # could also use desc(logFC)
  slice(1:5)
top_genes <- unique(c(top10$gene, bottom10$gene))
length(top_genes)
data_label <- diff_g1[diff_g1$gene %in% top_genes, ]
colors <- c("up" = "#fe8264", "down" = "#b5dbe6", "normal" = "#474747")
colors
# Inspect the data
head(diff_g1)
p <- ggplot(diff_g1, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(data = data_bg, shape = 21, fill = "#474747", color = "black", alpha = 0.8, size = 1.3, stroke = 0.3) +  # non-significant grey points (shape 21 = bordered circle; stroke = border width)
  geom_point(data = data_sig, aes(fill = color), shape = 21, color = "black", size = 2, stroke = 0.3) +  # significant points
  scale_fill_manual(values = colors,  # point colours
                    labels = c("up" = "up in MP", "down" = "up in P")  # custom legend labels
                    ) +
  geom_text_repel(data = data_label, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), size = 3,
                  force_pull = 0,         # label attraction set to 0 so labels are not pulled back to the points
                  point.padding = 0,      # minimum distance between a label and its point
                  box.padding = 1,        # minimum distance between a label and the points
                  min.segment.length = 0, # add a segment whenever its length is > 0
                  vjust = 0.5,            # right side of the label aligned to the position
                  segment.color = "grey20",
                  segment.size = 0.5,     # segment thickness
                  segment.alpha = 0.8,    # segment transparency
                  max.overlaps = Inf,
                  seed = 123, max.time = 1, max.iter = Inf,
                  nudge_y = 0.5, nudge_x = -0.3) +  # small label offset along y
  geom_hline(yintercept = -log10(0.05), color = "#474747", linewidth = 0.6, linetype = "dashed") +  # significance threshold line
  scale_y_continuous(expand = c(0, 0)) +  # points start right at the axis line
  labs(title = "DEGs of Tumor Core Macrophages: MP vs P", x = expression(log[2](FC)), y = expression(-log[10](adj.p))) +
  theme_classic() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
ggsave(filename = file.path(PROJECT_ROOT, "figures", "macrophage_DEG_tumor_core_MP_vs_P.pdf"), width = 8, height = 6, plot = p, bg = "white")
