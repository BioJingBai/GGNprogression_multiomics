# -----------------------------------------------------------------------------
# Canonical marker dot plot across the Level-4 cell types
#
# Purpose   : Compute the positive markers of each final_Level4_celltype and
#             draw the canonical-marker DotPlot (expression level + fraction of
#             expressing cells) used as the annotation reference panel.
#
# Inputs    : <PROJECT_ROOT>/processed/smi_annotated.qs
#
# Outputs   : <PROJECT_ROOT>/figures/canonical_marker_dotplot.pdf
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript 04_figureS13c_canonical_marker_dotplot.R
# -----------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Dependencies
## ---------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(qs)

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
dir.create(file.path(PROJECT_ROOT, "figures"), recursive = TRUE, showWarnings = FALSE)

## Load the annotated object produced by 00_setup_and_annotation.R
smi <- qread(file.path(PROJECT_ROOT, "processed", "smi_annotated.qs"))


## ---------------------------------------------------------------------------
## Marker reference
## ---------------------------------------------------------------------------
# EP: EPCAM ELF3
# Mix:
# Plasma: MZB1 JCHAIN IGHA1 IGHG1/2
# (macrophage core): LYZ CD68 C1QA C1QB C1QC APOE APOC1 CD163
# ma_SPP1: SPP1
# (endothelial core): VWF PECAM1 FLT1 PLVAP
# endothelia_Bro: MYC HSPG2 SPRY1 VWA1 HBEGF
# (T core): CD3D CD3E CD2 TRAC TRBC2 TIGIT
# CD4T: CD4
# NKT: NKG7 FGFBP2 IFNG
# mesothelial: KRT18 SLPI
# activated_B: CD69 CD83 EGR1
# endothelia_Lym: CCL21 PDPN
# ma_COX2: COX2
# NoFs: DCN IGFBP6 MFAP5
# SMC: ACTA2 MYL9 TAGLN
# CD8T: CD8A CD8B
# (B core): MS4A1 CD19 CD79A LTB
# Naive_B: TCL1A FCER2 IGHD
# ma_STAB1: STAB1
# DC: AIF1 CYBB CLEC9A LAMP3 CLEC4C CD207
# ma_CCL18: CCL18
# ma_MARCO: MARCO
# CAFs: TAGLN POSTN COL1A1 COL3A1 COL6A3 ACTA2 RGS5 CTHRC1


## ---------------------------------------------------------------------------
## Dot plot
## ---------------------------------------------------------------------------
# First look at the top markers
Idents(smi) <- "final_Level4_celltype"
DefaultAssay(smi) <- "RNA"
DEGs_allpos <- FindAllMarkers(smi, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
DEGs_allpos = DEGs_allpos[DEGs_allpos$p_val_adj < 0.05, ]
DEGs_allpos %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
all_genes <- c(
  # EP
  "EPCAM", "ELF3", "KRT18",
  "LYZ", "CD68", "C1QA", "C1QB", "C1QC", "APOE", "APOC1", "CD163",
  # ma_MARCO
  "MARCO",
  # ma_CCL18
  "CCL18",
  # ma_SPP1
  "SPP1",
  # ma_STAB1
  "STAB1",
  # ma_COX2
  "COX2",
  # DC
  "AIF1", "CYBB", "CLEC9A", "LAMP3", "CLEC4C", "CD207",
  "CD3D", "CD3E", "CD2", "TRAC", "TRBC2", "TIGIT",
  # CD4T
  "CD4",
  # CD8T
  "CD8A", "CD8B",
  # NKT
  "NKG7", "FGFBP2", "IFNG",
  "MS4A1", "CD19", "CD79A", "LTB",
  # activated_B
  "CD69", "CD83", "EGR1",
  # Naive_B
  "TCL1A", "FCER2", "IGHD",
  "VWF", "PECAM1", "FLT1", "PLVAP",
  # Plasma
  "MZB1", "JCHAIN", "IGKC", "IGHA1", "IGHG1", "IGHG2",
  # endothelia_Bro
  "MYC", "HSPG2", "SPRY1", "VWA1",
  # endothelia_Lym
  "CCL21", "PDPN",
  # CAFs
  "POSTN", "COL1A1", "COL3A1", "COL6A3", "ACTA2", "CTHRC1",
  # NoFs (normal fibroblasts)
  "DCN",
  # SMC
  "ACTA2", "MYL9", "TAGLN",
  # Mesothelial
  "SLPI"
)
all_genes <- unique(all_genes)
# Set the display order of the cell types
# Reverse the level order
level_order <- rev(c(
  "EP", "ma_MARCO", "ma_CCL18", "ma_SPP1", "ma_STAB1", "ma_COX2", "DC", "CD4T",
  "CD8T", "NKT", "activated_B", "Naive_B", "Plasma",
  "endothelia_Bro", "endothelia_Lym", "CAFs", "NoFs", "SMC", "mesothelial", "Mix"
))
smi@meta.data$final_Level4_celltype <- factor(
  smi@meta.data$final_Level4_celltype,
  levels = level_order
)
Idents(smi) <- "final_Level4_celltype"
p <- DotPlot(
  smi,
  features = all_genes,
  cols = c("#4575B4", "#F46D43"),
  dot.scale = 6
) +
  scale_color_gradientn(
    colours = c("#313695", "#74ADD1", "#FEE090", "#F46D43", "#A50026")
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 40, hjust = 1, vjust = 1, size = 12   # gene names tilted for readability
    ),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "white", colour = "white"),
    strip.text = element_text(size = 13, face = "bold")
  )
pdf(file = file.path(PROJECT_ROOT, "figures", "canonical_marker_dotplot.pdf"), width = 18, height = 6)
print(p)
dev.off()
