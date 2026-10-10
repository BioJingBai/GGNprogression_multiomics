# -----------------------------------------------------------------------------
# GSEA bar plots (GO-BP / KEGG / Hallmark)
#
# Purpose   : Select the top positive-score pathways per group from the enrichment result tables and draw one horizontal bar plot per gene-set collection.
#
# Inputs    : figures/GO_BP_plot.csv, figures/hallmark_plot.csv, figures/KEGG_plot.csv (written by script 06)
#
# Outputs   : figures/GSEA_GO.pdf, figures/GSEA_Hallmark.pdf, figures/GSEA_KEGG.pdf
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript figureS4bd_gsea_plots.R
# -----------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

GO_BP_plot<-read_csv(file.path(PROJECT_ROOT, "figures/GO_BP_plot.csv"))
GO_BP_plot_S<- dplyr::filter(GO_BP_plot,GO_BP_plot$comparison=="Svsother") |>
  dplyr::filter(score>0) |>
  arrange(desc(score))|>
  mutate(group="S")
GO_BP_plot_P<- dplyr::filter(GO_BP_plot,GO_BP_plot$comparison=="Pvsother")|>
  dplyr::filter(score>0) |>
  arrange(desc(score))|>
  mutate(group="P")
GO_BP_plot_MP<- dplyr::filter(GO_BP_plot,GO_BP_plot$comparison=="MPvsother")|>
  dplyr::filter(score>0) |>
  arrange(desc(score))|>
  mutate(group="MP")
GO_BP_plot_P<-GO_BP_plot_P[-c(3,4,5,7,8,9,14,15),]
GO_BP_plot_S<-GO_BP_plot_S |>head(3)
GO_BP_plot_P<-GO_BP_plot_P |>head(4)
GO_BP_plot_MP<-GO_BP_plot_MP |>head(4)
GO_BP_plot<-rbind(GO_BP_plot_S,GO_BP_plot_P,GO_BP_plot_MP)
colors <- c(
  "S"  = "#F3E6C8",  # cream apricot yellow (lighter than the previous yellow)
  "P"  = "#C7D2E1",  # macaron blue-grey
  "MP" = "#E6B8B7"   # macaron bean-paste pink
)
GO_BP_plot$group<-factor(GO_BP_plot$group,levels = c("MP","P","S"))
GO_BP_plot2 <- GO_BP_plot %>%
  arrange(group, desc(-log10(pvalue))) %>%    # sort within group by -logP (largest at the bottom)
  mutate(ID = factor(ID, levels = unique(ID)))
p <- ggplot(data = GO_BP_plot2, aes(x = -log10(pvalue), y = ID, fill = group)) +
  geom_bar(width = 0.7, stat = 'identity',color="black",) +
  geom_text(aes(x = 0.1, label = ID), size = 3.5, hjust = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(x = "-Log10(Pvalue)",
       y = "GO_BP Processing") +
  theme_classic() +
  theme(
    legend.position = 'none',
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1)
  )
ggsave(file.path(PROJECT_ROOT, "figures/GSEA_GO.pdf"), plot = p, width = 4.5, height = 3.5)

h_plot<-read_csv(file.path(PROJECT_ROOT, "figures/hallmark_plot.csv"))
h_plot_S<- dplyr::filter(h_plot,h_plot$comparison=="Svsother")|>
  dplyr::filter(score>0) |>
  arrange(desc(score))|>
  mutate(group="S")
h_plot_P<- dplyr::filter(h_plot,h_plot$comparison=="Pvsother")|>
  dplyr::filter(score>0) |>
  arrange(desc(score))|>
  mutate(group="P")
h_plot_MP<- dplyr::filter(h_plot,h_plot$comparison=="MPvsother")|>
  dplyr::filter(score>0) |>
  arrange(desc(score))|>
  mutate(group="MP")
h_plot_S<-h_plot_S |>head(3)
h_plot_P<-h_plot_P |>head(4)
h_plot_MP<-h_plot_MP |>head(4)
h_plot<-rbind(h_plot_S,h_plot_P,h_plot_MP)
colors <- c("S"  = "#F3E6C8", "P"  = "#C7D2E1", "MP" = "#E6B8B7")
h_plot$group<-factor(h_plot$group,levels = c("MP","P","S"))
h_plot2 <- h_plot %>%
  arrange(group, desc(-log10(pvalue))) %>%    # sort within group by -logP (largest at the bottom)
  mutate(ID = factor(ID, levels = unique(ID)))
p <- ggplot(data = h_plot2, aes(x = -log10(pvalue), y = ID, fill = group)) +
  geom_bar(width = 0.7, stat = 'identity',color="black",) +
  geom_text(aes(x = 0.1, label = ID), size = 3.5, hjust = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(x = "-Log10(Pvalue)",
       y = "Hallmark Pathway") +
  theme_classic() +
  theme(
    legend.position = 'none',
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1)
  )
ggsave(file.path(PROJECT_ROOT, "figures/GSEA_Hallmark.pdf"), plot = p, width = 4.5, height = 3.5)

KEGG_plot<-read_csv(file.path(PROJECT_ROOT, "figures/KEGG_plot.csv"))
KEGG_plot_S<- dplyr::filter(KEGG_plot,KEGG_plot$comparison=="Svsother")|>
  dplyr::filter(score>0) |>
  arrange(desc(score))|>
  mutate(group="S")
KEGG_plot_P<- dplyr::filter(KEGG_plot,KEGG_plot$comparison=="Pvsother")|>
  dplyr::filter(score>0) |>
  arrange(desc(score))|>
  mutate(group="P")
KEGG_plot_MP<- dplyr::filter(KEGG_plot,KEGG_plot$comparison=="MPvsother")|>
  dplyr::filter(score>0) |>
  arrange(desc(score))|>
  mutate(group="MP")
KEGG_plot_P<-KEGG_plot_P[-c(3,5,6,7,9,12),]
KEGG_plot_MP<-KEGG_plot_MP[-c(4,5,6,9,10),]
KEGG_plot_S<-KEGG_plot_S |>head(3)
KEGG_plot_P<-KEGG_plot_P |>head(4)
KEGG_plot_MP<-KEGG_plot_MP |>head(4)
KEGG_plot<-rbind(KEGG_plot_S,KEGG_plot_P,KEGG_plot_MP)
#colors <- c("S"="#F5E4BB","P"="#8491B4FF","MP"='#e15759')
colors <- c("S"  = "#F3E6C8", "P"  = "#C7D2E1", "MP" = "#E6B8B7")
KEGG_plot$group<-factor(KEGG_plot$group,levels = c("MP","P","S"))
KEGG_plot2 <- KEGG_plot %>%
  arrange(group, desc(-log10(pvalue))) %>%    # sort within group by -logP (largest at the bottom)
  mutate(ID = factor(ID, levels = unique(ID)))
p <- ggplot(data = KEGG_plot2, aes(x = -log10(pvalue), y = ID, fill = group)) +
  geom_bar(width = 0.7, stat = 'identity',color="black",) +
  geom_text(aes(x = 0.1, label = ID), size = 3.5, hjust = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(x = "-Log10(Pvalue)",
       y = "KEGG Processing") +
  theme_classic() +
  theme(
    legend.position = 'none',
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1)
  )
ggsave(file.path(PROJECT_ROOT, "figures/GSEA_KEGG.pdf"), plot = p, width = 4.5, height = 3.5)
