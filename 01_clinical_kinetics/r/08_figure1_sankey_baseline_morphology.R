# -----------------------------------------------------------------------------
# Sankey, baseline-composition and morphology panels
#
# Purpose   : Draw the first-to-last CT Sankey, the baseline nodule-type
#             composition by trajectory, and the morphology-feature panels
#             (including the standalone vascular-sign version).
#
# Inputs    : $PROJECT_ROOT/omics/supplementary_clinical_analysis/figure_assembly/rerun_999/input/TableS1_curated_binary_nodule_type_999.xlsx
#
# Outputs   : $PROJECT_ROOT/.../rerun_999/figures/*.pdf / *.png ;
#             $PROJECT_ROOT/.../rerun_999/ledger/panel_statistics.csv
#
# Usage     : LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript 08_figure1_sankey_baseline_morphology.R
# -----------------------------------------------------------------------------

if (!grepl("UTF-8", toupper(Sys.getlocale("LC_CTYPE")))) {
  cat("\n  LC_CTYPE is not UTF-8; the Chinese column headers of the source tables will not match. Run with:\n",
      "    LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript 08_figure1_sankey_baseline_morphology.R\n\n", sep = "")
}

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(patchwork); library(ggalluvial); library(scales)
})

.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}
HERE <- .script_dir()
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", "/path/to/GGN_progression_data")
PROJ <- sub("/$", "", PROJECT_ROOT)
FIG  <- file.path(PROJ, "omics/supplementary_clinical_analysis/figure_assembly/rerun_999/figures"); dir.create(FIG, showWarnings = FALSE)
RES  <- file.path(PROJ, "omics/supplementary_clinical_analysis/figure_assembly/rerun_999/ledger"); dir.create(RES, showWarnings = FALSE)
TABLES1 <- file.path(PROJ, "omics/supplementary_clinical_analysis/figure_assembly/rerun_999/input/TableS1_curated_binary_nodule_type_999.xlsx")

# ---- Colours: two palettes across the figure, no overlap --------------------
GRP_COL <- c("Stable" = "#4E79A7", "GGO-growth" = "#F28E2B", "SOLID-growth" = "#E15759")
GRP_LV  <- c("Stable", "GGO-growth", "SOLID-growth")
# Nodule type: only these two use this palette (light = ground glass, dark = part solid)
BIN_COL <- c("pGGN" = "#F7C3A4", "PSN" = "#AC7BA6")   # one step darker than the Sankey palette so the ribbons stay visible

theme_cc <- function(base = 13) {
  theme_classic(base_size = base) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = base + 1),
          axis.text.x   = element_text(size = base - 2, colour = "black"),
          axis.text.y   = element_text(size = base - 2, colour = "black"),
          axis.title    = element_text(size = base - 1, face = "bold"),
          legend.position = "none",
          axis.line     = element_line(linewidth = 0.6),
          axis.ticks    = element_line(linewidth = 0.6),
          plot.margin   = margin(6, 8, 6, 6))
}
save_panel <- function(p, name, w, h) {
  pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h); print(p); dev.off()
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  cat(sprintf("  [OK] %-40s %.1f x %.1f in\n", name, w, h))
}
pfmt <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))

STATS <- list()
rec <- function(panel, item, value)
  STATS[[length(STATS) + 1]] <<- data.frame(panel = panel, item = item,
                                            value = as.character(value), stringsAsFactors = FALSE)

# =============================================================================
# Read the data (same definition as the rebuild script: binary TableS1)
# =============================================================================
cat("\n================ read data ================\n")
t1 <- read_excel(TABLES1) %>% as.data.frame()
# Source column headers are in Chinese; map them to ASCII names once, here, so that
# every downstream reference is locale-independent (R cannot parse Chinese
# identifiers under a C locale).
names(t1)[names(t1) == "基线结节性质"] <- "baseline_appearance"   # baseline nodule type (0 = pGGN, 1 = PSN)
names(t1)[names(t1) == "末次结节性质"] <- "final_appearance"      # final-visit nodule type (0 = pGGN, 1 = PSN)
names(t1)[names(t1) == "血管穿行"]   <- "vascular_transit"      # vascular sign (0/1)
names(t1)[names(t1) == "贴近胸膜"]   <- "pleural_attachment"    # pleural attachment (0/1)
names(t1)[names(t1) == "边缘规则"]   <- "regular_margin"        # regular margin (0/1)
names(t1)[names(t1) == "分叶征"]     <- "lobulation"            # lobulation sign (0/1)
names(t1)[names(t1) == "毛刺征"]     <- "spiculation"           # spiculation sign (0/1)
names(t1)[names(t1) == "空泡征"]     <- "air_bubble"            # air-bubble sign (0/1)
names(t1)[names(t1) == "胸膜牵拉"]   <- "pleural_retraction"    # pleural retraction (0/1)
cat(sprintf("TableS1: %d rows x %d columns\n", nrow(t1), ncol(t1)))

t1 <- t1 %>%
  mutate(Group = factor(Group, levels = GRP_LV)) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Baseline_type = factor(ifelse(baseline_appearance == 0, "pGGN", "PSN"), levels = c("pGGN", "PSN")),
    Final_type    = factor(ifelse(final_appearance == 0, "pGGN", "PSN"), levels = c("pGGN", "PSN")),
    Vascular      = suppressWarnings(as.numeric(vascular_transit)),
    Pleural_tag   = suppressWarnings(as.numeric(pleural_attachment)),
    Margin_smooth = suppressWarnings(as.numeric(regular_margin)),
    Lobulation    = suppressWarnings(as.numeric(lobulation)),
    Spiculation   = suppressWarnings(as.numeric(spiculation)),
    Air_bubble    = suppressWarnings(as.numeric(air_bubble)),
    Pleural_pull  = suppressWarnings(as.numeric(pleural_retraction))
  ) %>%
  mutate(across(c(Vascular, Pleural_tag, Margin_smooth, Lobulation,
                  Spiculation, Air_bubble, Pleural_pull), ~ .x > 0,
                .names = "{.col}_pos"))

cat(sprintf("three groups: %s\n", paste(names(table(t1$Group)), table(t1$Group), sep = "=", collapse = " ")))
cat(sprintf("first-CT nodule type: %s\n",
            paste(names(table(t1$Baseline_type)), table(t1$Baseline_type), sep = "=", collapse = " ")))

# =============================================================================
# [d] Sankey - third axis recoloured blue/orange/red (first/second axes stay peach/violet)
# =============================================================================
cat("\n================ [d] Sankey ================\n")

df_d <- t1 %>% select(Baseline_type, Final_type, Group)
df_d <- to_lodes_form(df_d, axes = 1:3, id = "value")
df_d$stratum <- factor(df_d$stratum,
                       levels = c("pGGN", "PSN", "Stable", "GGO-growth", "SOLID-growth"))
stopifnot(!any(is.na(df_d$stratum)))

sankey_col <- c(BIN_COL, GRP_COL)     # pGGN/PSN peach+violet, the three trajectories blue/orange/red

p_d <- ggplot(df_d, aes(x = x, fill = stratum, label = stratum,
                        stratum = stratum, alluvium = value)) +
  geom_flow(width = 0.51, curve_type = "sine", alpha = 0.72, color = "white", linewidth = 0.08) +
  geom_stratum(width = 0.5, color = "grey35", linewidth = 0.25) +
  geom_text(stat = "stratum", size = 3.2, color = "black", lineheight = 0.88,
            aes(label = sprintf("%s\n%.1f%%", after_stat(stratum),
                                100 * after_stat(count) / nrow(t1)))) +
  scale_fill_manual(values = sankey_col) +
  scale_x_discrete(labels = c("First CT", "Last CT", "Trajectory"),
                   expand = expansion(add = c(0.08, 0.08))) +
  theme_void() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 13, face = "bold", colour = "black"),
        plot.margin = margin(6, 10, 6, 10))

save_panel(p_d, "Fig1d_sankey_first_to_last", 7.4, 5.6)
for (i in seq_len(nrow(df_d))) NULL
mig <- t1 %>% count(Baseline_type, Final_type, Group, name = "n") %>%
  mutate(pct_of_total = round(100 * n / sum(n), 1)) %>%
  arrange(Group, Baseline_type, Final_type)
print(as.data.frame(mig), row.names = FALSE)
for (i in seq_len(nrow(mig)))
  rec("d-sankey", paste(mig$Group[i], mig$Baseline_type[i], "→", mig$Final_type[i]),
      paste0(mig$n[i], " (", mig$pct_of_total[i], "% of 999)"))

# =============================================================================
# [e] baseline composition (flipped) - first-CT type -> trajectory
# =============================================================================
cat("\n================ [e] baseline composition (flipped) ================\n")

ct_e <- table(t1$Baseline_type, t1$Group)
chi_e <- suppressWarnings(chisq.test(ct_e))
cat(sprintf("  chi-square = %.2f, df = %d, p = %.3g\n",
            unname(chi_e$statistic), unname(chi_e$parameter), chi_e$p.value))
print(ct_e)
rec("e-baseline", "Chi-square p", formatC(chi_e$p.value, format = "e", digits = 3))

df_e <- as.data.frame(prop.table(ct_e, margin = 1)) %>%
  setNames(c("Type", "Group", "Prop")) %>%
  left_join(as.data.frame(ct_e) %>% setNames(c("Type", "Group", "n")), by = c("Type", "Group")) %>%
  mutate(Type  = factor(Type,  levels = c("pGGN", "PSN")),
         Group = factor(Group, levels = GRP_LV)) %>%
  arrange(Type, Group) %>%
  group_by(Type) %>% mutate(cumtop = cumsum(Prop)) %>% ungroup()

# The denominator (nodules of each first-CT type) goes into the x-axis labels
n_by_type <- df_e %>% group_by(Type) %>% summarise(N = sum(n), .groups = "drop")
xlabs <- sprintf("%s\n(n=%d)", n_by_type$Type, n_by_type$N)

p_e <- ggplot(df_e, aes(x = Type, y = Prop, fill = Group)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.5) +
  geom_text(aes(y = cumtop - Prop / 2, label = sprintf("%d (%.1f%%)", n, 100 * Prop)),
            size = 3.0, colour = "white", fontface = "bold") +
  scale_fill_manual(values = GRP_COL, name = NULL) +
  scale_x_discrete(labels = xlabs) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Trajectory outcome by nodule type at first CT", x = NULL, y = "Proportion",
       subtitle = sprintf("Chi-square test: p = %s", pfmt(chi_e$p.value))) +
  theme_cc(13) +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10.5, colour = "black"),
        plot.subtitle = element_text(hjust = 0.5, size = 10.5, colour = "grey25"))

save_panel(p_e, "Fig1e_nodule_type_first_CT", 5.2, 4.4)

for (tt in levels(df_e$Type)) {
  sub <- df_e %>% filter(Type == tt)
  rec("e-baseline", paste0("first-CT ", tt, " | n"), sum(sub$n))
  for (i in seq_len(nrow(sub)))
    rec("e-baseline", paste0("first-CT ", tt, " -> ", sub$Group[i]),
        sprintf("%d (%.1f%%)", sub$n[i], 100 * sub$Prop[i]))
}

# The two numbers used in the text: share of first-CT pGGN that grow / share of first-CT PSN that stay stable
prog_pggn <- df_e %>% filter(Type == "pGGN", Group != "Stable") %>% summarise(s = sum(Prop)) %>% pull(s)
stab_psn  <- df_e %>% filter(Type == "PSN",  Group == "Stable")  %>% summarise(s = sum(Prop)) %>% pull(s)
cat(sprintf("  * first-CT pGGN that grow: %.1f%% | first-CT PSN that stay stable: %.1f%%\n",
            100 * prog_pggn, 100 * stab_psn))
rec("e-baseline", "first-CT pGGN that grow (%)", sprintf("%.1f", 100 * prog_pggn))
rec("e-baseline", "first-CT PSN that stay stable (%)",     sprintf("%.1f", 100 * stab_psn))

# =============================================================================
# [S1f] morphology - recoloured to "one colour per group"
# =============================================================================
# Original: fill = Sign(No/Yes), blue = No / orange = Yes - which clashed with orange = GGO-growth.
# New version: stacking order Yes at the bottom, No at the top; both segments use
# **the group colour**, and the No segment fades to 18% opacity. So "bar colour =
# group" (one colour per group) and "solid height = positive rate", which stays
# readable. The legend lists only the three groups, one scale for the whole figure.
# =============================================================================
cat("\n================ [S1f] morphology ================\n")

MORPH <- list(
  list(var = "Vascular_pos",      lab = "Vascular sign"),
  list(var = "Pleural_tag_pos",   lab = "Pleural attachment"),
  list(var = "Margin_smooth_pos", lab = "Regular margin"),
  list(var = "Lobulation_pos",    lab = "Lobulation"),
  list(var = "Spiculation_pos",   lab = "Spiculation"),
  list(var = "Air_bubble_pos",    lab = "Air bubble"),
  list(var = "Pleural_pull_pos",  lab = "Pleural retraction"))

MORPH_PLOTS <- list()
for (k in seq_along(MORPH)) {
  v <- MORPH[[k]]$var; lab <- MORPH[[k]]$lab
  dd <- t1 %>% filter(!is.na(.data[[v]])) %>% mutate(Group = droplevels(Group))
  ct <- table(dd$Group, dd[[v]])
  ch <- suppressWarnings(chisq.test(ct))
  cat(sprintf("  %-20s n=%3d | chi2=%6.2f df=%d p=%s | %s\n",
              lab, nrow(dd), unname(ch$statistic), unname(ch$parameter),
              formatC(ch$p.value, format = "e", digits = 3),
              paste(rownames(ct), sprintf("%d/%d", ct[, "TRUE"], rowSums(ct)),
                    sep = "=", collapse = " ")))
  rec("S1f-morphology", paste0(lab, " | n per group (Present/Total)"),
      paste(rownames(ct), sprintf("%d/%d", ct[, "TRUE"], rowSums(ct)), sep = "=", collapse = "; "))
  rec("S1f-morphology", paste0(lab, " | Chi-square p"),
      formatC(ch$p.value, format = "e", digits = 3))

  dp <- dd %>%
    group_by(Group) %>%
    summarise(Present = sum(.data[[v]]), Total = n(), .groups = "drop") %>%
    mutate(Absent = Total - Present, Prop = Present / Total) %>%
    pivot_longer(c("Present", "Absent"), names_to = "Sign", values_to = "n") %>%
    mutate(Sign = factor(Sign, levels = c("Present", "Absent")),
           Group = factor(Group, levels = GRP_LV)) %>%
    arrange(Group, Sign) %>%
    group_by(Group) %>% mutate(cumtop = cumsum(n) / Total) %>% ungroup() %>%
    mutate(pct = 100 * n / Total)

  MORPH_PLOTS[[k]] <- ggplot(dp, aes(x = Group, y = n / Total, fill = Group, alpha = Sign)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.4) +
    # ggplot stacks the **first factor level at the top**, so Present (solid segment)
    # is at the top and Absent (faded) at the bottom. The label sits at the centre of
    # the solid segment = 1 - Prop/2; a segment below 7% cannot hold it, so the label
    # moves below the segment and the white text switches to dark (below it is the
    # faded segment, where white would disappear).
    geom_text(data = subset(dp, Sign == "Present"),
              aes(y = ifelse(pct < 7, 1 - n / Total - 0.012, 1 - (n / Total) / 2),
                  label = sprintf("%d(%.1f%%)", n, pct)),
              colour = ifelse(subset(dp, Sign == "Present")$pct < 7, "grey20", "white"),
              size = 2.4, fontface = "bold",
              vjust = ifelse(subset(dp, Sign == "Present")$pct < 7, 1, 0.5)) +
    scale_fill_manual(values = GRP_COL, name = NULL, guide = "none") +
    scale_alpha_manual(values = c("Present" = 1, "Absent" = 0.18),
                       name = NULL, guide = "none") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.06))) +
    scale_x_discrete(expand = expansion(add = 0.40)) +
    labs(title = lab, x = NULL, y = "Proportion",
         subtitle = sprintf("Chi-square test: p = %s", pfmt(ch$p.value))) +
    theme_cc(11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
          plot.title = element_text(size = 10.5),
          plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey25"))
}

p_morph <- wrap_plots(MORPH_PLOTS, ncol = 4) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") +
  theme(legend.position = "right")
save_panel(p_morph, "FigS1f_morphology_features", 13.2, 7.8)

# ---- S1i vascular sign: standalone version of cell a of the morphology panel ---
# Same data and definition (n=823, p=1.6e-05), rendered separately. Blue = No / orange = Yes.
p_s1i <- MORPH_PLOTS[[1]] +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
        plot.title    = element_text(size = 15),
        plot.subtitle = element_text(size = 11),
        plot.margin   = margin(8, 12, 8, 8))
save_panel(p_s1i, "FigS1i_vascular_sign", 4.2, 3.8)

# Narrow version: 2 columns x 4 rows - one tall strip in the main figure, and each cell is larger after scaling than in the 13.2-inch version
p_morph_c <- wrap_plots(MORPH_PLOTS, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") +
  theme(legend.position = "none")
save_panel(p_morph_c, "Fig1morph_compact", 7.4, 8.6)

# ---- Write out ----------------------------------------------------------
out <- do.call(rbind, STATS)
write.csv(out, file.path(RES, "panel_statistics.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("\nledger -> ledger/panel_statistics.csv  (%d rows)\n", nrow(out)))
