# -----------------------------------------------------------------------------
# Bulk SPP1 versus immune-checkpoint gene correlation
#
# Purpose   : Compute Spearman correlations between SPP1 expression and the exhaustion/checkpoint genes and draw the faceted scatter plot.
#
# Inputs    : objects from the upstream R session: mrna_count_tumor (created by script 03; this script transposes it in place). No files are read from disk.
#
# Outputs   : figures/SPP1_exhaustion_gene_correlation.pdf
#
# Usage     : PROJECT_ROOT=/path/to/data Rscript figureS20b_spp1_checkpoint_correlation.R
# -----------------------------------------------------------------------------

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")

mrna_count_tumor<-t(mrna_count_tumor)
exhaustion_genes <- c(
  "LAG3", "HAVCR2", "TNFRSF9",
  "TIGIT", "PDCD1"
)
genes_use <- c("SPP1", exhaustion_genes)
expr_mat <- mrna_count_tumor[, genes_use]
expr_mat<-as.data.frame(expr_mat)
cor_res <- map_df(
  exhaustion_genes,
  function(g) {
    ct <- cor.test(
      expr_mat$SPP1,
      expr_mat[[g]],
      method = "spearman",
      exact = FALSE
    )
    tibble(
      gene = g,
      rho = ct$estimate,
      p.value = ct$p.value,
      direction = ifelse(ct$estimate > 0, "Positive", "Negative")
    )
  }
)

# reshape expr_mat to long format for plotting
expr_long <- expr_mat %>%
  mutate(sample = rownames(expr_mat)) %>%  # add a sample column
  pivot_longer(
    cols = all_of(exhaustion_genes),
    names_to = "gene",
    values_to = "expr"
  )

# compute the Spearman correlation of each gene
cor_stats <- expr_long %>%
  group_by(gene) %>%
  summarise(
    rho = cor(expr, SPP1, method = "spearman"),
    p.value = cor.test(expr, SPP1, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    label = paste0("rho=", round(rho, 2), " ", signif)
  )

# plot
p <- ggplot(expr_long, aes(x = expr, y = SPP1)) +
  geom_point(size = 1.5, alpha = 0.7, color = "#2c7fb8") +
  geom_smooth(method = "lm", se = FALSE, color = "#f03b20", linetype = "dashed") +
  facet_wrap(~gene, scales = "free_x") +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  xlab("Gene expression") +
  ylab("SPP1 expression") +
  geom_text(
    data = cor_stats,
    aes(
      x = -Inf, y = Inf, 
      label = label
    ),
    hjust = -0.1, vjust = 1.2,
    inherit.aes = FALSE,
    size = 3.5
  )
pdf(file = file.path(PROJECT_ROOT, "figures/SPP1_exhaustion_gene_correlation.pdf"), width =6, height = 4)
print(p)
dev.off()
