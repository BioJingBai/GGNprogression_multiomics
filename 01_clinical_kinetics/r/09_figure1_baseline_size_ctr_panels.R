# -----------------------------------------------------------------------------
# Baseline size and CTR panels
#
# Purpose   : Draw the baseline CT panels in a consistent layout and style
#             (ggplot2 / theme_classic / scale_fill_tableau / ggpubr / ggalluvial
#             / theme_void, PDF output): Sankey, initial-status composition,
#             first growth event, whole-lesion diameter, solid-component volume,
#             CTR, vascular sign and morphology.
#
# Inputs    : $PROJECT_ROOT/omics/supplementary_clinical_analysis/figure_assembly/rerun_999/input/TableS1_curated_binary_nodule_type_999.xlsx
#             $PROJECT_ROOT/clinical_data/first_growth_events.xlsx  (panel e only)
#
# Outputs   : ./figures/<panel>.pdf (vector, for submission) + .png (300 dpi preview)
#             ./results/Figure1_panel_statistics.csv
#
# Usage     : LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript 09_figure1_baseline_size_ctr_panels.R
# -----------------------------------------------------------------------------

# ---- Encoding pre-check ------------------------------------------------------
# The source tables have Chinese column headers, so this script must run under a UTF-8 locale;
# otherwise R decodes the Chinese names as Latin-1 and read_excel fails with "zip file ... cannot be opened".
# Run: LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript 09_figure1_baseline_size_ctr_panels.R
if (!grepl("UTF-8", toupper(Sys.getlocale("LC_CTYPE")))) {
  cat("\n  current LC_CTYPE = ", Sys.getlocale("LC_CTYPE"), "\n",
      "    the Chinese column headers of the source tables will not match. Use instead:\n",
      "    LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript 09_figure1_baseline_size_ctr_panels.R\n\n", sep = "")
}

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)   # whole-cohort panel 1g: two stacked sections, box plot + zero-value strip
  library(ggpubr)
  library(ggalluvial)
  library(ggthemes)
  library(scales)
  library(stringr)
})

# ------------------------------------------------------- Paths (portable) -----
# Data paths resolve from PROJECT_ROOT; set it in the environment (see the header).
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

FIG  <- file.path(HERE, "figures"); dir.create(FIG,  showWarnings = FALSE)
RES  <- file.path(HERE, "results"); dir.create(RES,  showWarnings = FALSE)

TABLES1 <- file.path(PROJ, "omics/supplementary_clinical_analysis/figure_assembly/rerun_999/input/TableS1_curated_binary_nodule_type_999.xlsx")
CLINDAT <- file.path(PROJ, "clinical_data/first_growth_events.xlsx")

# ------------------------------------------------------------------ Global style
mycol <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
           "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF")

# Three-group palette: the first three ggthemes Tableau 10 colours (blue / orange /
# red), as used by scale_fill_tableau().
GRP_COL <- c("Stable" = "#4E79A7", "GGO-growth" = "#F28E2B", "SOLID-growth" = "#E15759")
GRP_LV  <- c("Stable", "GGO-growth", "SOLID-growth")

# Binary nodule type: solid orange / red
BIN_COL <- c("pGGN" = "#F28E2B", "PSN" = "#E15759")

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
  cat(sprintf("  [OK] %-34s %.1f x %.1f in\n", name, w, h))
}

# p-value formatting: scientific notation below 1e-3, otherwise 3 decimals, so that
# values 3 decimals cannot express (4.2e-4 -> "0.000") stay legible.
pfmt <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))

STATS <- list()
rec <- function(panel, item, value) {
  STATS[[length(STATS) + 1]] <<- data.frame(panel = panel, item = item, value = as.character(value),
                                            stringsAsFactors = FALSE)
}

KW <- function(df, y, panel, label) {
  f  <- as.formula(paste(y, "~ Group"))
  kt <- kruskal.test(f, data = df)
  ns <- df %>% group_by(Group) %>% summarise(n = sum(!is.na(.data[[y]])), .groups = "drop")
  cat(sprintf("  [%s] %s  KW chi2=%.2f df=%d p=%.3g | n: %s\n",
              panel, label, unname(kt$statistic), unname(kt$parameter), kt$p.value,
              paste(ns$Group, ns$n, sep = "=", collapse = " ")))
  rec(panel, paste0(label, " | n per group"), paste(ns$Group, ns$n, sep = "=", collapse = "; "))
  rec(panel, paste0(label, " | Kruskal-Wallis p"), formatC(kt$p.value, format = "e", digits = 3))
  kt
}

# =============================================================================
# Read the data
# =============================================================================
cat("\n================ read data ================\n")

t1 <- read_excel(TABLES1) %>% as.data.frame()
# Source column headers are in Chinese; map them to ASCII names once, here, so that
# every downstream reference is locale-independent (R cannot parse Chinese
# identifiers under a C locale).
names(t1)[names(t1) == "影像号"]                   <- "image_id"                # imaging accession number (the patient key)
names(t1)[names(t1) == "基线结节性质"]             <- "baseline_appearance"     # baseline nodule type (0 = pGGN, 1 = PSN)
names(t1)[names(t1) == "末次结节性质"]             <- "final_appearance"        # final-visit nodule type (0 = pGGN, 1 = PSN)
names(t1)[names(t1) == "VDT_2D优先3D补_天"]        <- "VDT_primary_days"        # primary VDT in days (2D first, 3D back-filled, growth-only)
names(t1)[names(t1) == "VDT_2D优先3D补_带符号_天"] <- "VDT_signed_days"         # the same with a sign; supplementary figures only
names(t1)[names(t1) == "Total_Volume_mm3_初次"]    <- "Total_Volume_mm3_first"  # whole-lesion volume at the first CT (mm3)
names(t1)[names(t1) == "Total_Volume_mm3_末次"]    <- "Total_Volume_mm3_last"   # whole-lesion volume at the last CT (mm3)
names(t1)[names(t1) == "Solid_Volume_mm3_初次"]    <- "Solid_Volume_mm3_first"  # solid-component volume at the first CT (mm3)
names(t1)[names(t1) == "Solid_Volume_mm3_末次"]    <- "Solid_Volume_mm3_last"   # solid-component volume at the last CT (mm3)
names(t1)[names(t1) == "CTR_3D_初次"]              <- "CTR_3D_first"            # solid/total volume ratio at the first CT
names(t1)[names(t1) == "CTR_3D_末次"]              <- "CTR_3D_last"             # solid/total volume ratio at the last CT
names(t1)[names(t1) == "血管穿行"]                 <- "vascular_transit"        # vascular sign (0/1)
names(t1)[names(t1) == "贴近胸膜"]                 <- "pleural_attachment"      # pleural attachment (0/1)
names(t1)[names(t1) == "边缘规则"]                 <- "regular_margin"          # regular margin (0/1)
names(t1)[names(t1) == "分叶征"]                   <- "lobulation"              # lobulation sign (0/1)
names(t1)[names(t1) == "毛刺征"]                   <- "spiculation"             # spiculation sign (0/1)
names(t1)[names(t1) == "空泡征"]                   <- "air_bubble"              # air-bubble sign (0/1)
names(t1)[names(t1) == "胸膜牵拉"]                 <- "pleural_retraction"      # pleural retraction (0/1)
names(t1)[names(t1) == "随访时长(天)"]             <- "follow_up_days"          # follow-up length from the first to the last CT (days)
names(t1)[names(t1) == "MDT_3D总质量_仅增长_天"]   <- "MDT_total_growing_days"  # mass doubling time of the whole lesion (days, growth-only)
names(t1)[names(t1) == "MDT_3D实性质量_仅增长_天"] <- "MDT_solid_growing_days"  # mass doubling time of the solid component (days, growth-only)
names(t1)[names(t1) == "年化CTR_3D"]               <- "annualized_CTR_3D"       # annualized 3D CTR change (/yr)
names(t1)[names(t1) == "3D来源"]                   <- "seg_source"              # which batch the 3D segmentation came from
cat(sprintf("TableS1: %d rows x %d columns | distinct image IDs %d | %s\n",
            nrow(t1), ncol(t1), length(unique(t1$image_id)),
            if ("Group" %in% names(t1)) paste(names(table(t1$Group)), table(t1$Group), sep = "=", collapse = " ") else "no Group column"))

# Group: unify to the full English name and fix the level order
t1 <- t1 %>%
  mutate(Group = factor(Group, levels = GRP_LV)) %>%
  filter(!is.na(Group))

# Binary nodule type: 0 = pGGN, 1 = PSN (the source table is already binarised)
t1 <- t1 %>%
  mutate(
    Baseline_type = factor(ifelse(baseline_appearance == 0, "pGGN", "PSN"), levels = c("pGGN", "PSN")),
    Final_type    = factor(ifelse(final_appearance == 0, "pGGN", "PSN"), levels = c("pGGN", "PSN")),
    VDT_main      = suppressWarnings(as.numeric(VDT_primary_days)),          # primary definition, growth-only
    VDT_signed    = suppressWarnings(as.numeric(VDT_signed_days)),    # signed; used only in the supplementary figures
    V0_3d         = suppressWarnings(as.numeric(Total_Volume_mm3_first)),
    V1_3d         = suppressWarnings(as.numeric(Total_Volume_mm3_last)),
    S0_3d         = suppressWarnings(as.numeric(Solid_Volume_mm3_first)),
    S1_3d         = suppressWarnings(as.numeric(Solid_Volume_mm3_last)),
    CTR0_3d       = suppressWarnings(as.numeric(CTR_3D_first)),
    CTR1_3d       = suppressWarnings(as.numeric(CTR_3D_last)),
    Vascular      = suppressWarnings(as.numeric(vascular_transit)),
# The remaining morphology features alongside the vascular sign (same imaging descriptions, same missingness source)
# NOTE: spiculation has one row with value 10 (should be 1), so every feature uses "> 0 = present" rather than == 1
    Pleural_tag   = suppressWarnings(as.numeric(pleural_attachment)),
    Margin_smooth = suppressWarnings(as.numeric(regular_margin)),
    Lobulation    = suppressWarnings(as.numeric(lobulation)),
    Spiculation   = suppressWarnings(as.numeric(spiculation)),
    Air_bubble    = suppressWarnings(as.numeric(air_bubble)),
    Pleural_pull  = suppressWarnings(as.numeric(pleural_retraction)),
    # The follow-up column name contains parentheses, so it must be back-quoted, or R parses it as a function call
    Years         = suppressWarnings(as.numeric(follow_up_days)) / 365.25
  )

# --- The analysis set shared by 1f / 1g / 1h ---
# SEG_OK: nodules with a 3D segmentation (the three panels share one analysis set, so they stay comparable)
t1 <- t1 %>%
  mutate(
    SEG_OK = !is.na(V0_3d),
    D0_3d  = ifelse(SEG_OK & V0_3d > 0, (6 * V0_3d / pi)^(1/3), NA_real_)
  )

# --- Data fix: 102/890 rows of the first-CT solid-volume column were divided by 1000 --------
# Symptom: those rows satisfy Solid_Volume / Total_Volume ~ 0.001 (exactly 1/1000),
#        while CTR_3D_first of the same row holds the "correct" value (the two differ by >= 5e-4).
#       Example: V0 = 9925.13, S0 = 9.9251 (= V0/1000), yet CTR = 1.0000.
#       The last-CT column has the same problem (160/889 rows).
# Diagnosis: the solid-volume column carries a unit error on those rows (mm3 -> cm3); the CTR and total-volume columns are reliable.
# Fix: back-derive the solid volume as CTR x total volume, replacing only the rows where the stored column disagrees with CTR;
#       rows that agree with CTR keep their stored value (CTR was stored rounded, and back-deriving would add rounding error).
# Effect: the solid-component diameter of panel 1g and the CTR of panel 1h are both based on the corrected quantity.
t1 <- t1 %>%
  mutate(
    S0_ok  = !is.na(S0_3d) & !is.na(CTR0_3d) & V0_3d > 0 &
             abs(S0_3d / V0_3d - CTR0_3d) < 5e-4,
    S0_rep = ifelse(S0_ok, S0_3d, CTR0_3d * V0_3d),
# Equivalent spherical diameter of the solid component (0 when the solid volume is 0)
    SD0_rep = ifelse(!SEG_OK | is.na(S0_rep), NA_real_,
                     ifelse(S0_rep > 0, (6 * S0_rep / pi)^(1/3), S0_rep))
  )

cat(sprintf("solid-volume fix: %d/%d rows disagreed with CTR and were replaced with CTR x total volume\n",
            sum(!t1$S0_ok & !is.na(t1$S0_3d) & !is.na(t1$CTR0_3d) & t1$V0_3d > 0),
            sum(!is.na(t1$S0_3d) & !is.na(t1$CTR0_3d) & t1$V0_3d > 0)))
rec("1g/1h", "Solid_Volume rows replaced with CTR x total volume",
    sum(!t1$S0_ok & !is.na(t1$S0_3d) & !is.na(t1$CTR0_3d) & t1$V0_3d > 0))

cat(sprintf("1f/1g/1h shared analysis set (3D segmentation available): n = %d | pGGN %d / PSN %d\n",
            sum(t1$SEG_OK), sum(t1$SEG_OK & t1$Baseline_type == "pGGN"),
            sum(t1$SEG_OK & t1$Baseline_type == "PSN")))

cat(sprintf("plotted: %d nodules | %s\n", nrow(t1),
            paste(names(table(t1$Group)), table(t1$Group), sep = "=", collapse = " ")))

# =============================================================================
# Fig 1b - Sankey: first CT -> last CT -> group (binary pGGN / PSN)
#   to_lodes_form + geom_flow(curve_type = "sine") + theme_void
# =============================================================================
cat("\n================ Fig 1b  Sankey ================\n")

df_b <- t1 %>% select(Baseline_type, Final_type, Group)
key  <- c("Baseline_type", "Final_type", "Group")
df_b <- to_lodes_form(df_b, axes = 1:length(key), id = "value")
# stratum is the per-axis value and may only take these 5 levels (writing transition combinations turns them all NA)
df_b$stratum <- factor(df_b$stratum,
                       levels = c("pGGN", "PSN", "Stable", "GGO-growth", "SOLID-growth"))

# Sankey palette: pastel, distinct from the group colours and the binary-type colours
col_b <- c("pGGN" = "#FFD7C5", "PSN" = "#C496B8",
           "Stable" = "#FFF4D1", "GGO-growth" = "#F1DCD4", "SOLID-growth" = "#DCE8E2")

p_b <- ggplot(df_b, aes(x = x, fill = stratum, label = stratum,
                        stratum = stratum, alluvium = value)) +
  geom_flow(width = 0.51, curve_type = "sine", alpha = 0.6, color = "white", linewidth = 0.08) +
  geom_stratum(width = 0.5, color = "grey35", linewidth = 0.25) +
# "name + percentage" on two lines inside each block. All three axes contain every
# nodule, so the denominator is the total nodule count.
  geom_text(stat = "stratum", size = 3.2, color = "black", lineheight = 0.88,
            aes(label = sprintf("%s\n%.1f%%", after_stat(stratum),
                                100 * after_stat(count) / nrow(t1)))) +
  scale_fill_manual(values = col_b) +
  scale_x_discrete(labels = c("First CT", "Last CT", "Group"),
                   expand = expansion(add = c(0.08, 0.08))) +
  theme_void() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 13, face = "bold", colour = "black"),
        plot.margin = margin(6, 10, 6, 10))

save_panel(p_b, "Fig1b_Sankey_first_last_group", 7.4, 5.6)

# Transition counts (for the caption and the text cross-check)
mig <- t1 %>% count(Baseline_type, Final_type, Group, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% arrange(Group, Baseline_type, Final_type)
write.csv(mig, file.path(RES, "Fig1b_transition_counts.csv"), row.names = FALSE)
cat("  transition counts:\n"); print(as.data.frame(mig), row.names = FALSE)
for (i in seq_len(nrow(mig)))
  rec("1b", paste(mig$Group[i], mig$Baseline_type[i], "→", mig$Final_type[i]), paste0(mig$n[i], " (", mig$pct[i], "%)"))

# =============================================================================
# Fig 1c - five continuous-kinetics panels (VDT / MDT total / MDT solid / annualized CTR / follow-up years)
#   violin + white box overlay (log10 panels), jittered points (zero-inflated /
#   long-tailed panels), KW p in the top-left, pairwise Mann-Whitney + BH stars,
#   x-axis labels (n = entering the analysis / nodules in that group).
#   Data: the primary-definition columns of TableS1, growth-only nodules.
# =============================================================================
cat("\n================ Fig 1c  five continuous-kinetics panels ================\n")

GROUPS <- GRP_LV
# Pairwise Mann-Whitney with BH correction (within each panel's 3 comparisons; rank tests are invariant to the log transform)
pair_bh <- function(df, y, panel, label) {
  d <- df %>% filter(!is.na(.data[[y]]))
  # Denominator = all nodules of that trajectory (including those missing this metric), not the number entering the test
  tot <- df %>% count(Group, name = "total")
  prs <- list(c(1, 2), c(1, 3), c(2, 3))
  pv  <- sapply(prs, function(k) {
    sub <- d %>% filter(Group %in% GROUPS[k])
    suppressWarnings(wilcox.test(sub[[y]] ~ droplevels(sub$Group))$p.value)
  })
  padj <- p.adjust(pv, method = "BH")
  ns <- d %>% count(Group, name = "n")
  kt <- kruskal.test(as.formula(paste(y, "~ Group")), data = d)
  cat(sprintf("  [%s] %s | n: %s | KW p=%.3g | BH p: %s\n", panel, label,
              paste(ns$Group, ns$n, sep = "=", collapse = " "), kt$p.value,
              paste(sprintf("%.3g", padj), collapse = " / ")))
  rec(panel, paste0(label, " | n per group"), paste(ns$Group, ns$n, sep = "=", collapse = "; "))
  rec(panel, paste0(label, " | KW p"), formatC(kt$p.value, format = "e", digits = 3))
  rec(panel, paste0(label, " | pairwise BH p (S-GP, S-AP, GP-AP)"),
      paste(sprintf("%.3g", padj), collapse = "; "))
  list(d = d, ns = ns, tot = tot, padj = padj, kw = kt$p.value)
}

star_of <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

# Draw one continuous-kinetics panel
kinetic_panel <- function(df, y, ylab, panel_id, letter, use_violin = TRUE) {
  st <- pair_bh(df, y, panel_id, ylab)
  d  <- st$d
  nlab <- setNames(sprintf("%s\n(n=%d/%d)", st$ns$Group, st$ns$n,
                           st$tot$total[match(st$ns$Group, st$tot$Group)]), st$ns$Group)

# Bracket positions are always computed in plotting coordinates first and converted back to data units:
#   on a log10 panel the plotting coordinate is a log value while annotate() takes a data value, so a 10^ back-transform is needed;
#   on a linear panel the two are identical.
  if (use_violin) {
    # Violin + white box, log10 y axis (VDT/MDT span 3-4 orders of magnitude; a linear axis flattens the three groups)
    lim <- range(d[[y]])
    br  <- 10^(floor(log10(lim[1])):ceiling(log10(lim[2])))
    amp <- diff(log10(lim)) * 0.10              # vertical spacing between brackets (in log10 units)
    ys  <- log10(lim[2]) + amp * (0.6 + 0:2)    # log10 coordinates
    ys_d <- 10^ys                               # back to data units
    h_d  <- 10^(ys + 0.25 * amp) - ys_d         # bracket leg height (data units)
    p <- ggplot(d, aes(x = Group, y = .data[[y]], fill = Group)) +
      geom_violin(width = 0.65, alpha = 0.5, colour = "black", linewidth = 0.6, trim = TRUE) +
      geom_boxplot(width = 0.14, fill = "white", alpha = 0.9, colour = "black",
                   linewidth = 0.7, outlier.shape = NA) +
      scale_y_log10(breaks = br, labels = scales::comma_format())
  } else {
    # Zero-inflated / long-tailed metrics: jittered points + white box (instead of a violin KDE)
    set.seed(20260912)
    rg  <- range(d[[y]], na.rm = TRUE)
    amp <- diff(rg) * 0.06
    ys_d <- rg[2] * 0.86 + amp * (0:2)
    h_d  <- amp * 0.25
    p <- ggplot(d, aes(x = Group, y = .data[[y]], fill = Group)) +
      geom_jitter(width = 0.13, size = 0.35, alpha = 0.35, colour = "grey25") +
      geom_boxplot(width = 0.20, fill = "white", alpha = 0.9, colour = "black",
                   linewidth = 0.7, outlier.shape = NA) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "grey50") +
      scale_y_continuous()
  }

  # Pairwise BH stars (brackets share the data scale)
  seg <- do.call(rbind, lapply(1:3, function(i) {
    k <- list(c(1, 2), c(1, 3), c(2, 3))[[i]]
    data.frame(x1 = k[1], x2 = k[2], y = ys_d[i], lab = star_of(st$padj[i]))
  }))

  p +
    scale_fill_manual(values = GRP_COL) +
    scale_x_discrete(labels = nlab) +
    annotate("segment", x = seg$x1, xend = seg$x1, y = seg$y, yend = seg$y + h_d,
             linewidth = 0.5, colour = "black") +
    annotate("segment", x = seg$x2, xend = seg$x2, y = seg$y, yend = seg$y + h_d,
             linewidth = 0.5, colour = "black") +
    annotate("segment", x = seg$x1, xend = seg$x2, y = seg$y + h_d,
             yend = seg$y + h_d, linewidth = 0.5, colour = "black") +
    annotate("text", x = (seg$x1 + seg$x2) / 2, y = seg$y + h_d * 1.25,
             label = seg$lab, size = 2.7) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.4, size = 2.6,
             colour = "dimgrey", label = sprintf("KW p = %.2e", st$kw)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.35, vjust = 3.6, size = 5,
             fontface = "bold", label = letter) +
    labs(x = NULL, y = ylab) +
    theme_cc(9.5) +
    theme(axis.text.x = element_text(size = 7.0, lineheight = 0.95),
          axis.title.y = element_text(size = 8.5, face = "plain"),
          plot.margin = margin(18, 10, 4, 6))
}

df_kin <- t1 %>% filter(!is.na(Group))
panel_c <- list(
  kinetic_panel(df_kin, "VDT_main",                     "VDT (days)",              "1c-A", "A", TRUE),
  kinetic_panel(df_kin, "MDT_total_growing_days",         "MDT total (days)",        "1c-B", "B", TRUE),
  kinetic_panel(df_kin, "MDT_solid_growing_days",       "MDT solid (days)",        "1c-C", "C", TRUE),
  kinetic_panel(df_kin, "annualized_CTR_3D",            "Annualized CTR change (/yr)", "1c-D", "D", FALSE),
  kinetic_panel(df_kin, "Years",                        "Follow-up (years)",       "1c-E", "E", FALSE)
)

# The metric columns must be numeric: converted once on df_kin
# Note: in patchwork 1.3.0 plot_annotation(theme=) raises "annotation$theme is not a valid theme",
# and `& theme()` does not work on patchwork objects in this version; the caption style comes from the global theme instead.
# Every panel carries its own theme_cc(), so this theme_set does not affect them.
old_theme <- theme_get()
theme_set(theme_get() + theme(plot.caption = element_text(size = 5.4, colour = "grey35", hjust = 0)))
p_c <- patchwork::wrap_plots(panel_c, ncol = 3) +
  patchwork::plot_annotation(
    caption = paste(
      "Kruskal-Wallis across the three trajectories; pairwise two-sided Mann-Whitney U with Benjamini-Hochberg",
      "adjustment across the three comparisons within each panel; *** p<0.001, ** p<0.01, * p<0.05, ns not significant.",
      "VDT and mass doubling times are shown on a log10 axis and include growing nodules only; all statistics computed on complete data.",
      "n = nodules entering the analysis / nodules in that trajectory.", sep = "\n"))

save_panel(p_c, "Fig1c_continuous_kinetics", 8.4, 5.4)

# =============================================================================
# Fig 1d - first-CT nodule-type composition (binary)
#   stacked proportion bars + chi-square p
# =============================================================================
cat("\n================ Fig 1d  Initial status ================\n")

ct_d <- table(t1$Group, t1$Baseline_type)
chi_d <- suppressWarnings(chisq.test(ct_d))
cat(sprintf("  chi-square = %.2f, df = %d, p = %.3g\n",
            unname(chi_d$statistic), unname(chi_d$parameter), chi_d$p.value))
print(ct_d)
rec("1d", "Chi-square p", formatC(chi_d$p.value, format = "e", digits = 3))
rec("1d", "n per group", paste(rownames(ct_d), rowSums(ct_d), sep = "=", collapse = "; "))

df_d <- as.data.frame(prop.table(ct_d, margin = 1)) %>%
  setNames(c("Group", "Type", "Prop")) %>%
  left_join(as.data.frame(ct_d) %>% setNames(c("Group", "Type", "n")), by = c("Group", "Type")) %>%
  mutate(Group = factor(Group, levels = GRP_LV))

p_d <- ggplot(df_d, aes(x = Group, y = Prop, fill = Type)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.4) +
  # The count label is sized for the longest string "395(68.7%)": about 47 pt at size 3.0, wider than the
  # 44 pt available in a 0.72-wide block, so its left end was clipped by the axis. Reduced to size 2.4 and the panel widened, keeping the whole string inside the block.
  geom_text(aes(label = sprintf("%d(%.1f%%)", n, 100 * Prop)),
            position = position_stack(vjust = 0.5), size = 2.4, colour = "white", fontface = "bold") +
  scale_fill_manual(values = BIN_COL, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "The initial status of GGNs", x = NULL, y = "Proportion",
       subtitle = sprintf("Chi-square test: p = %s", pfmt(chi_d$p.value))) +
  theme_cc(13) +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
        plot.subtitle = element_text(hjust = 0.5, size = 10.5, colour = "grey25"))

save_panel(p_d, "Fig1d_initial_status", 4.9, 4.2)

# =============================================================================
# Fig 1e - timing and magnitude of the first growth event
#   Data limitation: the first growth event is recorded only for nodules observed to grow, so n < the whole cohort
# =============================================================================
cat("\n================ Fig 1e  First growth event ================\n")

cd <- read_excel(CLINDAT, sheet = "Sheet1") %>% as.data.frame()
# The clinical table also has Chinese column headers; map them once, here.
names(cd)[names(cd) == "第一增长间隔"] <- "first_growth_interval"   # months from the first presentation to the first growth event
names(cd)[names(cd) == "大小增长"]     <- "size_change"            # size change at that first growth event (mm)
cat(sprintf("clinical data: %d rows | group != NG: %d\n", nrow(cd), sum(cd$group != "NG")))

df_e <- cd %>%
  filter(group != "NG") %>%
  transmute(ID,
            Group = recode(group, "S" = "Stable", "P" = "GGO-growth", "MP" = "SOLID-growth"),
            first_interval = suppressWarnings(as.numeric(first_growth_interval)),
            size_change    = suppressWarnings(as.numeric(size_change))) %>%
  filter(!is.na(first_interval), !is.na(size_change)) %>%
  mutate(Group = factor(Group, levels = GRP_LV))

cat(sprintf("  usable points n = %d | %s\n", nrow(df_e),
            paste(names(table(df_e$Group)), table(df_e$Group), sep = "=", collapse = " ")))
print(df_e %>% group_by(Group) %>%
        summarise(n = n(), median_months = median(first_interval),
                  median_mm = median(size_change), .groups = "drop") %>% as.data.frame(), row.names = FALSE)
rec("1e", "n per group (with first-growth event)",
    paste(names(table(df_e$Group)), table(df_e$Group), sep = "=", collapse = "; "))

p_e <- ggplot(df_e, aes(x = first_interval, y = size_change)) +
  geom_hline(yintercept = 2, color = "gray50", linetype = "dashed", linewidth = 0.6) +
  geom_point(aes(fill = Group), shape = 21, color = "black", size = 3, stroke = 0.5) +
  scale_fill_manual(values = GRP_COL, name = NULL) +
  labs(title = "Timing and magnitude of the first growth event",
       x = "Months from the first presentation",
       y = "The changes in size (mm)") +
  theme_minimal(base_size = 12) +
  theme(axis.text   = element_text(face = "bold", family = "sans", colour = "black"),
        axis.title  = element_text(face = "bold", family = "sans", size = 12),
        legend.position = "top",
        legend.text = element_text(size = 10),
        plot.title  = element_text(hjust = 0.5, face = "bold", size = 12.5),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA))

save_panel(p_e, "Fig1e_first_growth_event", 4.6, 4.4)

# =============================================================================
# Fig 1f-1i - baseline CT imaging features (whole cohort, 3D quantitative definition)
#     f equivalent spherical diameter of the whole lesion
#     g solid-component volume (PSN, log axis)
#     h CTR_3D (whole cohort)
#     i vascular sign
#   CTR is the direct ratio of the two 3D volumes, so solid <= total holds by
#   construction and the y axis never exceeds 1.
# =============================================================================
cat("\n================ Fig 1f–1i  Baseline imaging features ================\n")

box_panel <- function(df, y, title, ylab, name, panel_id, w = 3.1, h = 3.5,
                      logy = FALSE, zero_note = FALSE, zero_strip = FALSE,
                      set_tag = NULL) {
  d  <- df %>% filter(!is.na(.data[[y]]))
# all statistics are computed on the complete data

  # Log panels: nodules with zero solid volume cannot sit on a log axis (43% of the Stable group).
  # They are removed from the **display** only; the KW and pairwise tests still use the complete data, and the number excluded goes in the subtitle.
  dp     <- if (logy) d %>% filter(.data[[y]] > 0) else d
  n_excl <- nrow(d) - nrow(dp)
  if (n_excl > 0)
    rec(panel_id, paste0(title, " | not shown on the log axis (zero, still counted in the KW and pairwise tests)"),
        sprintf("%d/%d (%.1f%%)", n_excl, nrow(d), 100 * n_excl / nrow(d)))

  # The widest pair (1,3) goes on top - the usual way to nest brackets.
  CMP <- list(c(1, 2), c(2, 3), c(1, 3))

  # Pairwise Wilcoxon on the complete data, annotated with the uncorrected p:
  # ggpubr's label = "p.signif" reads the p column rather than p.adj. The
  # Holm-corrected values are written to the results table as well.
  pv <- vapply(CMP, function(k)
    suppressWarnings(wilcox.test(d[[y]][d$Group == GRP_LV[k[1]]],
                                 d[[y]][d$Group == GRP_LV[k[2]]])$p.value), numeric(1))
  pv_adj <- p.adjust(pv, method = "holm")
  star <- function(p) ifelse(p <= 1e-4, "****", ifelse(p <= 1e-3, "***",
                       ifelse(p <= 1e-2, "**",  ifelse(p <= 0.05, "*", "ns"))))
  rec(panel_id, paste0(title, " | pairwise Wilcoxon raw p (1v2 / 2v3 / 1v3)"),
      paste(formatC(pv, format = "e", digits = 3), collapse = "; "))
  rec(panel_id, paste0(title, " | the same, Holm-corrected (not used for the stars)"),
      paste(formatC(pv_adj, format = "e", digits = 3), collapse = "; "))

  # How high the brackets hang: above the **drawn upper whisker**, not above the data maximum.
  # With outlier.shape = NA the points beyond 1.5*IQR are not drawn at all, so positioning by max leaves the brackets
  # hanging far above the whisker (clearest in 1h: the whisker top is 0.52 but the brackets start at 1.07, all blank in between,
  # and the y axis is stretched to 1.3 - this is the "brackets hang too high" report).
  # On a log panel the box statistics are computed in log10 space (ggplot2 transforms before it stats),
  # so the whisker top is also found in log10 space and then back-transformed with 10^.
  tf <- if (logy) function(v) log10(v) else function(v) v
  bk <- if (logy) function(v) 10^v      else function(v) v
  whisk_top <- function(v) {
    t <- tf(v); q <- quantile(t, c(.25, .75), na.rm = TRUE)
    min(max(t, na.rm = TRUE), q[2] + 1.5 * diff(q))
  }
  sp <- split(dp[[y]], droplevels(dp$Group))
  sp <- sp[vapply(sp, length, integer(1)) > 0]
  ttop  <- max(vapply(sp, whisk_top, numeric(1)))
  tbase <- tf(min(dp[[y]], na.rm = TRUE))
  # About 0.16 in per level, enough for a 9-pt star at size 3.2. The two-section panel squeezes the box area to
  # 2.45 in, leaving about 9 pt per level with the stars touching the brackets, so the spacing is one notch wider.
  tstep <- (ttop - tbase) * (if (zero_strip) 0.145 else 0.11)
  tpos  <- ttop + tstep * seq_along(CMP)
  seg <- data.frame(
    x     = vapply(CMP, `[`, numeric(1), 1),
    xend  = vapply(CMP, `[`, numeric(1), 2),
    ybar  = bk(tpos),
    yfoot = bk(tpos - tstep * 0.20),
    ytxt  = bk(tpos + tstep * 0.08),
    lab   = star(pv)
  )

  # Log panels must account for the excluded zeros. The text has to be short - the panel is only 2.5 in wide,
  # the subtitle is drawn centred, and ggplot neither wraps nor clips, so anything too wide is cut mid-character
  # ("102/343 lesions at 0 mm3 not shown" left only 2 px of its right edge).
  sub <- sprintf("Kruskal-Wallis p = %s", pfmt(kt$p.value))
  # When sibling panels share a title (whole-cohort vs PSN solid volume), this line distinguishes the analysis sets
  if (!is.null(set_tag)) sub <- sprintf("%s\n%s", sub, set_tag)
  if (n_excl > 0 && !zero_strip)
    sub <- sprintf("%s\nlog scale; %d/%d at 0 mm³ not shown",
                   sub, n_excl, nrow(d))
  if (n_excl > 0 && zero_strip)
    # The lower strip already draws the zero part explicitly, so this only needs to say which part the box covers
    sub <- sprintf("%s\nlog scale; boxes, %d/%d at volume > 0",
                   sub, nrow(dp), nrow(d))

  # CTR panel: CTR is defined for **every** nodule (pGGN is 0 by definition) and the zero values are part of the
  # result, so they cannot be dropped with a log axis as with solid volume. The zero shares 67% / 44% / 21%
  # are exactly the gradient and go in the subtitle, so a reader of the boxes alone does not think all three groups sit at 0.
  if (zero_note) {
    zc <- vapply(GRP_LV, function(g) sum(d[[y]][d$Group == g] == 0, na.rm = TRUE), numeric(1))
    zn <- vapply(GRP_LV, function(g) sum(!is.na(d[[y]][d$Group == g])),           numeric(1))
    sub <- sprintf("%s\nCTR = 0: %s", sub, paste0(zc, "/", zn, collapse = " · "))
    rec(panel_id, paste0(title, " | CTR = 0 per group (n / total) (",
                         paste(GRP_LV, collapse = " / "), ")"),
        paste0(zc, "/", zn, collapse = "; "))
  }

  p <- ggplot(dp, aes(x = Group, y = .data[[y]], fill = Group)) +
    geom_boxplot(width = 0.6, alpha = 0.85, outlier.shape = NA, linewidth = 0.4) +
    scale_fill_manual(values = GRP_COL) +
    geom_segment(data = seg, inherit.aes = FALSE,
                 aes(x = x,     xend = x,     y = yfoot, yend = ybar), linewidth = 0.4) +
    geom_segment(data = seg, inherit.aes = FALSE,
                 aes(x = xend,  xend = xend,  y = yfoot, yend = ybar), linewidth = 0.4) +
    geom_segment(data = seg, inherit.aes = FALSE,
                 aes(x = x,     xend = xend,  y = ybar,  yend = ybar), linewidth = 0.4) +
    geom_text(data = seg, inherit.aes = FALSE, vjust = 0, size = 3.2,
              aes(x = (x + xend) / 2, y = ytxt, label = lab)) +
    labs(title = title, x = NULL, y = ylab, subtitle = sub) +
    theme_cc(11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
          plot.subtitle = element_text(hjust = 0.5, size = 8.0, colour = "grey25",
                                       lineheight = 1.1))

  # The upper limit is derived from the bracket heights (4.15 levels).
  # NOTE: this must never be written as scale limits = : ggplot drops out-of-range points **before** computing the box statistics,
  #    measured on 1h, three CTR points exceed 0.76, and with limits set the SOLID-growth Q3 moved from
  #    0.2240 to 0.2107 and the upper whisker from 0.5100 to 0.4600 - the box itself was changed.
  #    coord_cartesian only zooms the panel, leaving the data alone, and the box statistics are bit-identical (measured).
  ytop <- bk(ttop + tstep * 4.15)          # just above the highest star
  ybot <- bk(tbase - tstep * 0.45)         # a little room below so the box does not touch the axis
  sc <- if (logy) {
    brk <- 10^(seq(floor(log10(min(dp[[y]]))), ceiling(log10(max(dp[[y]])))))
    brk <- brk[brk >= min(dp[[y]]) & brk <= max(dp[[y]])]
    mnr <- as.vector(outer(2:9, brk))
    mnr <- mnr[mnr >= min(dp[[y]]) & mnr <= max(dp[[y]])]
    scale_y_log10(breaks = brk, minor_breaks = mnr,
                  labels = scales::label_number(big.mark = ","))
  } else {
    scale_y_continuous()
  }
  p <- p + sc + coord_cartesian(ylim = c(ybot, ytop))

  # --- Zero-value strip (whole-cohort definition) ---
  # Solid volume is exactly 0 in most nodules (503/890 = 56.5% of the whole cohort) and cannot be drawn on a log axis.
  # Dropping the zeros would leave a figure representing only 43.5% of the cohort, yet "is there any solid component" is exactly
  # where the three groups differ most. So a second section is stacked below the boxes showing = 0 and > 0;
  # all 890 nodules of the whole cohort appear, the upper half answering "how much when solid is present"
  # and the lower half "how many have none at all" - the full distribution needs both.
  if (zero_strip) {
    # The strip shows only the share with any solid component rather than a 100% stacked bar: the panel is under 1.6 in wide,
    # each bar is only 0.35 in wide while "367 (67%)" already takes 0.5 in at 8 pt, so forcing it into the block
    # bleeds out both sides (the first version cut the left half of "367"). The fix is a single bar with the label above it.
    ds <- d %>%
      mutate(has = .data[[y]] > 0) %>%
      group_by(Group) %>%
      summarise(n = sum(has), N = n(), Prop = n / N, .groups = "drop") %>%
      mutate(Group = factor(Group, levels = GRP_LV))
    rec(panel_id, paste0(title, " | solid volume = 0 per group (n / total) (",
                         paste(GRP_LV, collapse = " / "), ")"),
        paste(sprintf("%d/%d (%.1f%%)", ds$N - ds$n, ds$N, 100 * (1 - ds$Prop)),
              collapse = "; "))

    ps <- ggplot(ds, aes(x = Group, y = Prop, fill = Group)) +
      geom_col(width = 0.62, alpha = 0.85, colour = "white", linewidth = 0.4) +
      geom_text(aes(label = sprintf("%d (%.0f%%)", n, 100 * Prop)),
                vjust = -0.45, size = 2.3, colour = "grey20") +
      scale_fill_manual(values = GRP_COL, guide = "none") +
      scale_y_continuous(labels = percent_format(accuracy = 1),
                         # the strip is only 1 in tall, so five tick labels would crowd; cut to three
                         breaks = c(0, 0.5, 1), limits = c(0, 1),
                         expand = expansion(mult = c(0, 0.20))) +
      labs(x = NULL, y = NULL, subtitle = "Any solid component on CT") +
      theme_cc(11) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
            plot.subtitle = element_text(hjust = 0.5, size = 8.2, colour = "grey25"))

    # The upper section gives up its x labels so the strip carries them, which makes the two discrete axes align exactly
    p <- (p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                    plot.title = element_text(size = 11.5))) /
         ps + plot_layout(heights = c(2.45, 1))
  }

  save_panel(p, name, w, h)
  invisible(p)
}

# --- 1f equivalent spherical diameter of the whole lesion (all nodules with a 3D segmentation)
box_panel(t1 %>% filter(SEG_OK), "D0_3d", "GGN diameter", "Diameter (mm)",
          "Fig1f_GGN_diameter", "1f")

# --- 1g solid-component volume (restricted to baseline PSN: a log axis cannot hold 0) ----------
# A pure ground-glass nodule has zero solid volume by definition and cannot sit on a log axis, so 1g is restricted to baseline PSN -
# only there is the solid component a measurable continuous quantity. The number excluded is given in the subtitle.
t1_psn <- t1 %>% filter(SEG_OK, Baseline_type == "PSN")
cat(sprintf("\n  1g analysis set = baseline PSN with a 3D segmentation: n = %d | %s\n", nrow(t1_psn),
            paste(names(table(t1_psn$Group)), table(t1_psn$Group), sep = "=", collapse = " ")))
rec("1g", "analysis set (baseline PSN with 3D segmentation)", nrow(t1_psn))

# --- 1g solid-component volume, reported as a volume (no conversion to an
#     equivalent spherical diameter). On a linear axis the three boxes collapse into
#     a line near y = 0 (medians 1.3 / 23.5 / 91.1 mm3 against a maximum of 6659),
#     so only a log axis shows the gradient.
box_panel(t1_psn, "S0_rep", "Solid component volume",
          expression("Volume ("*mm^3*")"), "Fig1g_solid_volume", "1g", logy = TRUE)

# --- 1g-full solid-component volume, whole-cohort version (a second output; it does
#     not overwrite the panel above). Solid volume is defined for every nodule, as
#     CTR is, so the panel is not restricted to baseline PSN. Whole cohort n = 890,
#     of which 503 (56.5%) are exactly 0 and cannot sit on a log axis, so the panel
#     uses the two-section "box plot (> 0) + zero strip (all 890)" form, see
#     box_panel(zero_strip = TRUE). KW and pairwise tests run on all 890.
box_panel(t1 %>% filter(SEG_OK), "S0_rep", "Solid component volume",
          expression("Volume ("*mm^3*")"), "Fig1g_full_solid_volume", "1g-full",
          logy = TRUE, set_tag = "all lesions with 3D segmentation")

# --- 1h CTR (3D) - whole-cohort definition ---
# CTR is defined for every nodule (pGGN is 0 by definition), so 1h uses every nodule
# with a 3D segmentation (n = 890) without the PSN restriction of 1g. Group medians
# 0 / 0.006 / 0.060, Q3 0.006 / 0.040 / 0.189, zero shares 67% / 44% / 21%. CTR is
# the direct ratio of the 3D volumes, so solid <= total and the y axis never exceeds 1.
box_panel(t1 %>% filter(SEG_OK), "CTR0_3d", "CTR", "CTR (solid / total)",
          "Fig1h_CTR", "1h", zero_note = TRUE)

cat("\n  [1h] CTR out-of-range check (must be 0): ")
v <- t1$CTR0_3d[t1$SEG_OK]
cat(sprintf("n=%d, min=%.4f, max=%.4f,  count > 1 = %d\n",
            sum(!is.na(v)), min(v, na.rm = TRUE), max(v, na.rm = TRUE), sum(v > 1, na.rm = TRUE)))
rec("1h", "CTR out-of-range count (>1, must be 0)", sum(v > 1, na.rm = TRUE))

# --- 1i vascular sign
cat("\n  [1i] Vascular sign\n")
df_i <- t1 %>% filter(!is.na(Vascular)) %>%
  mutate(Vascular = factor(ifelse(Vascular == 1, "Yes", "No"), levels = c("No", "Yes")))
ct_i <- table(df_i$Group, df_i$Vascular)
chi_i <- suppressWarnings(chisq.test(ct_i))
cat(sprintf("  chi-square = %.2f, df = %d, p = %.3g\n",
            unname(chi_i$statistic), unname(chi_i$parameter), chi_i$p.value))
print(ct_i)
rec("1i", "Chi-square p", formatC(chi_i$p.value, format = "e", digits = 3))
rec("1i", "n per group", paste(rownames(ct_i), rowSums(ct_i), sep = "=", collapse = "; "))

df_i2 <- as.data.frame(prop.table(ct_i, margin = 1)) %>%
  setNames(c("Group", "Sign", "Prop")) %>%
  left_join(as.data.frame(ct_i) %>% setNames(c("Group", "Sign", "n")), by = c("Group", "Sign")) %>%
  mutate(Group = factor(Group, levels = GRP_LV))

p_i <- ggplot(df_i2, aes(x = Group, y = Prop, fill = Sign)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.4) +
  # The count label must stay inside the block: white text that leaves the block lands on white and disappears.
  # The label is wide, so the font is smaller and the panel wider, keeping the whole string inside the block.
  geom_text(aes(label = sprintf("%d(%.1f%%)", n, 100 * Prop),
                colour = Sign),   # dark text on light, white text on dark
            position = position_stack(vjust = 0.5), size = 2.4, fontface = "bold") +
  scale_colour_manual(values = c("No" = "white", "Yes" = "#7a3d05"), guide = "none") +
  # Original Fig 1j colours: No = blue, Yes = orange (the first two Tableau 10 colours), unchanged
  scale_fill_manual(values = c("No" = "#4E79A7", "Yes" = "#F28E2B"), name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Vascular sign on CT", x = NULL, y = "Proportion",
       subtitle = sprintf("Chi-square test: p = %s", pfmt(chi_i$p.value))) +
  theme_cc(11) +
  theme(legend.position = "right", legend.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey25"))

save_panel(p_i, "Fig1i_vascular_sign", 4.2, 3.8)

# =============================================================================
# 1k - morphology signs (7 panels)
#   1i draws only the vascular sign; the other 6 signs of TableS1 are drawn here in the
#   same layout (100% stacked bars by group + chi-square p), combined into panels a-g.
#   The only difference from 1i: positive counts sit above the bars rather than inside - lobulation / pleural retraction /
#   spiculation have positive rates as low as 1.4%-16.9%, so the label does not fit; all 7 panels move it up for consistency.
#   NOTE: these signs share the missingness of the vascular sign (~120-130/999), and it differs by group (see the ledger),
#   so each panel carries a "V5 batch only" sensitivity p for review.
# =============================================================================
cat("\n================ Fig 1k  Morphology features ================\n")

MORPH <- list(
  list(var = "Vascular",      lab = "Vascular sign on CT"),
  list(var = "Pleural_tag",   lab = "Pleural attachment"),
  list(var = "Margin_smooth", lab = "Regular margin"),
  list(var = "Lobulation",    lab = "Lobulation"),
  list(var = "Spiculation",   lab = "Spiculation"),
  list(var = "Air_bubble",    lab = "Air bubble"),
  list(var = "Pleural_pull",  lab = "Pleural retraction"))

MORPH_PLOTS <- list()
for (k in seq_along(MORPH)) {
  v <- MORPH[[k]]$var; lab <- MORPH[[k]]$lab
  # Presence is always "> 0": spiculation has one row with value 10 (should be 1), which == 1 would miss
  dd <- t1 %>% filter(!is.na(.data[[v]])) %>%
    # Labels are No / Yes.
    mutate(Sign = factor(ifelse(.data[[v]] > 0, "Yes", "No"),
                         levels = c("No", "Yes")),
           Group = droplevels(Group))
  ct <- table(dd$Group, dd$Sign)
  ch <- suppressWarnings(chisq.test(ct))
  cat(sprintf("  [1k] %-20s n=%3d | chi2=%6.2f df=%d p=%s | %s\n",
              lab, nrow(dd), unname(ch$statistic), unname(ch$parameter),
              formatC(ch$p.value, format = "e", digits = 3),
              paste(rownames(ct), sprintf("%d/%d", ct[, "Yes"], rowSums(ct)),
                    sep = "=", collapse = " ")))
  rec("1k", paste0(lab, " | n per group (Present/Total)"),
      paste(rownames(ct), sprintf("%d/%d", ct[, "Yes"], rowSums(ct)), sep = "=", collapse = "; "))
  rec("1k", paste0(lab, " | Chi-square p"), formatC(ch$p.value, format = "e", digits = 3))
  rec("1k", paste0(lab, " | missing count"), sprintf("%d/%d", sum(is.na(t1[[v]])), nrow(t1)))

  # Sensitivity: keep only the V5 batch (where these signs are missing in just 2.2%) to see whether the gradient survives
  dd5 <- dd %>% filter(seg_source == "V5(训练+内部测试)")
  if (nlevels(droplevels(dd5$Group)) == 3 && all(table(droplevels(dd5$Group), dd5$Sign) > 0)) {
    ch5 <- suppressWarnings(chisq.test(table(droplevels(dd5$Group), dd5$Sign)))
    rec("1k", paste0(lab, " | sensitivity: chi-square p, V5 batch only"),
        sprintf("%s (n=%d)", formatC(ch5$p.value, format = "e", digits = 3), nrow(dd5)))
  }

  dp <- as.data.frame(prop.table(ct, margin = 1)) %>%
    setNames(c("Group", "Sign", "Prop")) %>%
    left_join(as.data.frame(ct) %>% setNames(c("Group", "Sign", "n")),
              by = c("Group", "Sign")) %>%
    mutate(Group = factor(Group, levels = GRP_LV),
           Sign  = factor(Sign,  levels = c("No", "Yes"))) %>%
    # Cumulative top of the stacked bar: No at the bottom, Yes on top, so cumsum follows the Sign factor order
    arrange(Group, Sign) %>%
    group_by(Group) %>%
    mutate(cumtop = cumsum(Prop)) %>%
    ungroup()

  # Same style as the single vascular-sign panel: palette (No blue / Yes orange),
  # legend name, counts inside the blocks, the y-axis "Proportion" and the
  # "Chi-square test: p = ..." subtitle, arranged 4x2 on one page.
  MORPH_PLOTS[[k]] <- ggplot(dp, aes(x = Group, y = Prop, fill = Sign)) +
    geom_col(width = 0.80, colour = "white", linewidth = 0.4) +
    # White text that leaves a block lands on white and disappears, so light blocks use dark text and dark blocks white
    # A block below 7% cannot hold the label: a character is about 0.06 in tall while 1% is only 0.019 in,
    # so centring it overflows both ways - the top half turns white on the blue block and vanishes, the bottom half is cut by the panel edge
    # (this is how lobulation at 1.4% and pleural retraction at 2.7% were cut). The overflowing ones move above the block instead.
    geom_text(data = subset(dp, Sign == "Yes"),
              aes(y = ifelse(Prop < 0.07, cumtop + 0.012, cumtop - Prop / 2),
                  label = sprintf("%d(%.1f%%)", n, 100 * Prop), colour = Sign),
              size = 2.4, fontface = "bold", vjust = ifelse(
                subset(dp, Sign == "Yes")$Prop < 0.07, 0, 0.5)) +
    geom_text(data = subset(dp, Sign == "No"),
              aes(y = cumtop - Prop / 2, label = sprintf("%d(%.1f%%)", n, 100 * Prop),
                  colour = Sign),
              size = 2.4, fontface = "bold") +
    scale_colour_manual(values = c("No" = "white", "Yes" = "#7a3d05"), guide = "none") +
    scale_fill_manual(values = c("No" = "#4E79A7", "Yes" = "#F28E2B"), name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.06))) +
    # Discrete-axis padding tightened to 0.40: in a 4-column layout each panel is
    # about 2.7 in wide and the default 0.6 padding squeezes the bars below 0.5 in,
    # where "448(87.3%)" does not fit
    scale_x_discrete(expand = expansion(add = 0.40)) +
    labs(title = lab, x = NULL, y = "Proportion",
         subtitle = sprintf("Chi-square test: p = %s", pfmt(ch$p.value))) +
    theme_cc(11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
          plot.title = element_text(size = 10.5),
          plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey25"))
}

# The `&` operator (batch-applying a theme) does not work in this patchwork version,
# so `+ theme()` is used instead; it applies to the whole assembled figure.
p_morph <- wrap_plots(MORPH_PLOTS, ncol = 4) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") +
  theme(legend.position = "right", legend.text = element_text(size = 10))
save_panel(p_morph, "Fig1k_morphology_features", 13.2, 7.8)

# =============================================================================
# Summary output
# =============================================================================
cat("\n================ summary ================\n")

stat_df <- do.call(rbind, STATS)
write.csv(stat_df, file.path(RES, "Figure1_panel_statistics.csv"), row.names = FALSE)
cat(sprintf("statistics -> results/Figure1_panel_statistics.csv  (%d rows)\n", nrow(stat_df)))

# Per-panel analysis-set export (for review)
write.csv(t1 %>% select(image_id, Group, Baseline_type, Final_type, VDT_main, VDT_signed,
                        V0_3d, S0_3d, S0_ok, S0_rep, CTR0_3d, CTR1_3d,
                        D0_3d, SD0_rep, Vascular, Years),
          file.path(RES, "Fig1_analysis_set_full_cohort.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(df_e, file.path(RES, "Fig1e_analysis_set.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("\nfigures written:\n"); print(sort(list.files(FIG, pattern = "\\.pdf$")))
cat("\ndone.\n")
