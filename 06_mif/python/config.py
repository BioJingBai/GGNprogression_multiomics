#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF pipeline - shared configuration
#
# Purpose   : Input and output paths, the marker definitions of the two antibody
#             panels (P6 and P7), group order, labels and colours, figure defaults
#             and the Kruskal-Wallis plus pairwise statistical annotation helper.
#
# Inputs    : None. Every path is derived from PROJECT_ROOT:
#
# Outputs   : None.
#
# Usage     : import config   (or: from config import *)
#             Set PROJECT_ROOT to the directory holding the HALO exports before
#             running the pipeline.
# -----------------------------------------------------------------------------

import os

# ========================
# Paths
# ========================
# Root that contains: the per-panel HALO single-cell coordinate exports,
# the sample-to-group mapping workbook, and the working directory
# "mIF_Validation_Results" (intermediate tables + figures).
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", "/path/to/GGN_progression_data")

P6_BASE = os.path.join(PROJECT_ROOT, "mif_batch1_cells")
P7_BASE = os.path.join(PROJECT_ROOT, "mif_batch2_cells")
MAPPING_FILE = os.path.join(PROJECT_ROOT, "sample_id_map.xlsx")

# All outputs (figures, tables, stats) are written under OUTPUT_DIR.
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "mIF_Validation_Results")
OUTPUT_FIGURES = os.path.join(OUTPUT_DIR, "figures")
OUTPUT_DATA = os.path.join(OUTPUT_DIR, "data")
OUTPUT_STATS = os.path.join(OUTPUT_DIR, "stats")

# Vendor (company) reference statistics files — READ-ONLY inputs.
VENDOR_DIR = os.path.join(OUTPUT_DIR, "vendor_reference")
VENDOR_P6_STATS = os.path.join(VENDOR_DIR, "vendor_cell_stats_batch1.xlsx")
VENDOR_P7_STATS = os.path.join(VENDOR_DIR, "vendor_cell_stats_batch2.xlsx")

# ========================
# Groups and colors
# ========================
GROUP_ORDER = ["Stable", "Growth", "Solid"]
GROUP_LABELS = {"Stable": "Stable", "Growth": "GGO-growth", "Solid": "SOLID-growth"}
GROUP_COLORS = {"Stable": "#8DA0CB", "Growth": "#66C2A5", "Solid": "#FC8D62"}
GROUP_COLORS_LIGHT = {"Stable": "#C6DBEF", "Growth": "#B2E2D2", "Solid": "#FDD0A2"}

# ========================
# P6 Panel markers
# ========================
P6_MARKER_COLS = [
    "CD68 Positive Classification", "PD1 Positive Classification",
    "CD20 Positive Classification", "SPP1 Positive Classification",
    "GZMB Positive Classification", "CD8 Positive Classification",
    "CD4 Positive Classification", "PANCK Positive Classification",
    "aSMA Positive Classification"
]

P6_COMPOSITE_COLS = [
    "SPP1+CD68+", "PD1+CD8+", "PD1+CD4+", "PD1+CD4+CD8+",
    "GZMB+CD4+", "GZMB+CD8+", "GZMB+CD4+CD8+"
]

P6_ALL_MARKER_COLS = P6_MARKER_COLS + P6_COMPOSITE_COLS

P6_MARKER_NAMES = {
    "CD68 Positive Classification": "CD68+ Mac",
    "PD1 Positive Classification": "PD1+",
    "CD20 Positive Classification": "CD20+ B",
    "SPP1 Positive Classification": "SPP1+",
    "GZMB Positive Classification": "GZMB+",
    "CD8 Positive Classification": "CD8+ T",
    "CD4 Positive Classification": "CD4+ T",
    "PANCK Positive Classification": "PANCK+",
    "aSMA Positive Classification": "aSMA+",
    "SPP1+CD68+": "SPP1+CD68+ Mac",
    "PD1+CD8+": "PD1+CD8+ T",
    "PD1+CD4+": "PD1+CD4+ T",
    "PD1+CD4+CD8+": "PD1+CD4+CD8+",
    "GZMB+CD4+": "GZMB+CD4+ T",
    "GZMB+CD8+": "GZMB+CD8+ T",
    "GZMB+CD4+CD8+": "GZMB+CD4+CD8+",
}

# ========================
# P7 Panel markers
# ========================
P7_MARKER_COLS = [
    "CD68 Positive Classification", "CD31 Positive Classification",
    "CA9 Positive Classification", "CD20 Positive Classification",
    "CCR7 Positive Classification", "PDL1 Positive Classification",
    "CCL19 Positive Classification", "HLA-DR Positive Classification",
    "PANCK Positive Classification", "CD3 Positive Classification"
]

P7_COMPOSITE_COLS = [
    "CA9+CD31+", "PDL1+CD31+", "CA9+PDL1+CD31+", "HLA-DR+CD68+"
]

P7_ALL_MARKER_COLS = P7_MARKER_COLS + P7_COMPOSITE_COLS

P7_MARKER_NAMES = {
    "CD68 Positive Classification": "CD68+ Mac",
    "CD31 Positive Classification": "CD31+ Endo",
    "CA9 Positive Classification": "CA9+ Hypoxia",
    "CD20 Positive Classification": "CD20+ B",
    "CCR7 Positive Classification": "CCR7+",
    "PDL1 Positive Classification": "PDL1+",
    "CCL19 Positive Classification": "CCL19+",
    "HLA-DR Positive Classification": "HLA-DR+",
    "PANCK Positive Classification": "PANCK+",
    "CD3 Positive Classification": "CD3+ T",
    "CA9+CD31+": "CA9+CD31+ HypVes",
    "PDL1+CD31+": "PDL1+CD31+",
    "CA9+PDL1+CD31+": "CA9+PDL1+CD31+",
    "HLA-DR+CD68+": "HLA-DR+CD68+ Mac",
}

# ========================
# Coords columns
# ========================
COORD_COLS = ["XMin", "XMax", "YMin", "YMax"]
LABEL_COL = "Classifier Label"

# ========================
# Figure settings
# ========================
FIG_DPI = 300
FIG_FORMAT = "pdf"
FIG_SIZE_SINGLE = (5, 5)
FIG_SIZE_WIDE = (8, 5)
FIG_SIZE_SQUARE = (6, 6)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': FIG_DPI,
    'savefig.dpi': FIG_DPI,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

# Statistical annotation function
def get_stat_annotation(data_dict, group_order):
    """Run Kruskal-Wallis + pairwise Mann-Whitney (BH-FDR), return annotation dict."""
    from scipy.stats import kruskal, mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    import numpy as np

    groups = {g: data_dict[g] for g in group_order if g in data_dict and len(data_dict[g]) > 0}
    if len(groups) < 2:
        return {}

    values_list = [v for v in groups.values()]
    kw_stat, kw_p = kruskal(*values_list)

    pairs = []
    pvals = []
    pair_labels = []
    gkeys = list(groups.keys())
    for i in range(len(gkeys)):
        for j in range(i+1, len(gkeys)):
            g1, g2 = gkeys[i], gkeys[j]
            try:
                _, p = mannwhitneyu(groups[g1], groups[g2], alternative='two-sided')
                pairs.append((g1, g2))
                pvals.append(p)
                pair_labels.append((GROUP_LABELS.get(g1, g1), GROUP_LABELS.get(g2, g2)))
            except:
                pass

    if len(pvals) > 0:
        reject, pvals_corr, _, _ = multipletests(pvals, method='fdr_bh')
    else:
        pvals_corr = pvals

    annotations = {}
    for (g1, g2), (l1, l2), p in zip(pairs, pair_labels, pvals_corr):
        if p < 0.001:
            sig = "***"
        elif p < 0.01:
            sig = "**"
        elif p < 0.05:
            sig = "*"
        else:
            sig = "ns"
        annotations[(g1, g2)] = {"p": p, "sig": sig}

    return {"kw_p": kw_p, "pairwise": annotations}

print("Config loaded: P6 markers=%d, P7 markers=%d" % (len(P6_ALL_MARKER_COLS), len(P7_ALL_MARKER_COLS)))
