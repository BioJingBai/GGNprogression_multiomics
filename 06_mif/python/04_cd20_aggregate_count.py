#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF step 4 - CD20+ aggregate count per section
#
# Purpose   : Boxplot of the vendor-scored CD20+ aggregate count per section, with
#             the sections that have no aggregate scored 0 so that all 32 cohort
#             sections are shown.
#
# Inputs    : <VENDOR_DIR>/vendor_cd20_aggregate_counts_batch2.xlsx
#             <PROJECT_ROOT>/sample_id_map.xlsx, sheet 'batch2_map'
#             <OUTPUT_DATA>/P7_TC_metrics.csv   (cohort section list)
#
# Outputs   : <FIGS_DIR>/09_CD20_Aggregates/Fig9_Aggregate_Count.{pdf,png}
#
# Usage     : PROJECT_ROOT=/path/to/data python 04_cd20_aggregate_count.py
# -----------------------------------------------------------------------------

import os, re, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import figstyle as fs
warnings.filterwarnings('ignore')

np.random.seed(42)  # reproducibility (scatter jitter)

GROUP_ORDER = ['Stable', 'Growth', 'Solid']

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 11,
    'axes.titlesize': 14, 'axes.labelsize': 12, 'axes.titleweight': 'bold',
    'axes.labelweight': 'bold', 'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

PROJECT_ROOT = os.environ.get("PROJECT_ROOT", "/path/to/GGN_progression_data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "mIF_Validation_Results")
DATA_DIR = os.path.join(OUTPUT_DIR, "data")
FIGS_DIR = os.path.join(OUTPUT_DIR, "figures")
VENDOR_DIR = os.path.join(OUTPUT_DIR, "vendor_reference")
MAPPING_FILE = os.path.join(PROJECT_ROOT, "sample_id_map.xlsx")

print("Loading metrics ...")
print("Data loaded.")


def save_fig(fig, section, figname):
    d = os.path.join(FIGS_DIR, section)
    os.makedirs(d, exist_ok=True)
    fig.savefig(os.path.join(d, f"{figname}.pdf"), dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(d, f"{figname}.png"), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def get_group_data(df, col):
    dd = {}
    for g in GROUP_ORDER:
        vals = df[df['group'] == g][col].dropna().values
        if len(vals) >= 2:
            dd[g] = vals
    return dd


def r_boxplot(ax, df, col, title, ylabel="Percentage (%)"):
    """Publication-style tight boxplot panel: delegates drawing and the
    annotation layout to figstyle.panel.

    The statistical convention and the values are unchanged — grouping is still
    decided by get_group_data() (dropna + n >= 2 per group) and the test is
    still Kruskal-Wallis + two-sided Mann-Whitney U + Bonferroni correction;
    only the drawing style and the panel size change.
    """
    dd = get_group_data(df, col)
    present = [g for g in GROUP_ORDER if g in dd]
    if len(present) < 2:
        return None
    fs.panel(ax, dd, present, style='box', title=title, ylabel=ylabel)
    return True


# ============================================================
# CD20+ aggregates per section
# ============================================================
print("\n" + "=" * 40 + "\nCD20+ aggregates per section")
agg_file = os.path.join(VENDOR_DIR, "vendor_cd20_aggregate_counts_batch2.xlsx")
agg_raw = pd.read_excel(agg_file)
agg_raw.columns = ['ImageTag', 'Count']


# Two vendor file names carry typos that also occur in this xlsx:
# BT007050 -> true ID BT0007050; S0419341 -> true ID S0419314.
ID_ALIASES = {'BT007050': 'BT0007050', 'S0419341': 'S0419314'}


def extract_base_id(tag):
    """'...-20260716-BT0017501-5-MERGE.ome.tiff' -> 'BT0017501' (alias-normalized)."""
    m = re.search(r'-((?:BT|TS|S)\d+)(?:-\d+)?-M', str(tag), re.I)
    if not m:
        return None
    b = m.group(1)
    return ID_ALIASES.get(b, b)


# Fig9 works on the P7 panel (CD20+ aggregates): use the P7 mapping sheet,
# which includes a SOLID section that has P7 data only (vendor-confirmed).
map_df = pd.read_excel(MAPPING_FILE, sheet_name='batch2_map')
map_df['Group_clean'] = map_df['Group'].str.strip().replace({'Satble': 'Stable'})
# 'image_tag' holds the bare ID except one section that keeps its '-7'
# suffix; strip any trailing '-<digits>' so both sides compare as bare base IDs.
map_df['id_base'] = map_df['image_tag'].astype(str).str.strip().str.replace(r'-\d+$', '', regex=True)
id2g = dict(zip(map_df['id_base'], map_df['Group_clean']))

agg_raw['base_id'] = agg_raw['ImageTag'].apply(extract_base_id)
agg_raw['group'] = agg_raw['base_id'].map(id2g)
agg_all = agg_raw.dropna(subset=['group']).copy()

# Rows present in the vendor xlsx but absent from the mapping are dropped and
# reported explicitly.
dropped = agg_raw[agg_raw['group'].isna()]
print(f"  Aggregate rows mapped: {len(agg_all)}; explicitly excluded (not in mapping): {len(dropped)}")
for tag in dropped['ImageTag'].tolist():
    print(f"    - {tag}")
print("  Groups: " + ", ".join(f"{g}={int((agg_all['group'] == g).sum())}" for g in GROUP_ORDER))

# Sections with no circled aggregate are scored 0 (standard for a count metric).
# The section list is the cohort list of the P7 panel table, so every section is
# shown.
cohort_ids = pd.read_csv(os.path.join(DATA_DIR, 'P7_TC_metrics.csv'))['sample_id'].astype(str).str.replace(r'-\d+$', '', regex=True)
_have = set(agg_all['base_id'])
zero_added = []
for sid in cohort_ids:
    if sid in _have:
        continue
    g = id2g.get(sid)
    if g is None:
        print(f"    WARNING: {sid} has P7 data but is not in the mapping — group not invented, skipped")
        continue
    agg_all = pd.concat([agg_all, pd.DataFrame([{'base_id': sid, 'group': g, 'Count': 0}])],
                        ignore_index=True)
    zero_added.append((sid, g))
print(f"  Sections without aggregates included as 0: {len(zero_added)}")
for sid, g in zero_added:
    print(f"    + {sid} ({g}) -> Count = 0")
print("  Groups: " + ", ".join(f"{g}={int((agg_all['group'] == g).sum())}" for g in GROUP_ORDER))

# r_boxplot keeps every point visible as jitter and annotates the KW and
# Bonferroni pairwise p-values on the figure.
fig, ax = plt.subplots(figsize=fs.FIG_PANEL, facecolor='white')
r_boxplot(ax, agg_all, 'Count', 'CD20+ Aggregates', 'Count')
# Figure note: the zero-fill rule.
ax.annotate('No-aggregate sections scored 0; n = 32',
            xy=(0, -0.22), xycoords='axes fraction', ha='left', va='top',
            fontsize=5.8, color='0.25')
save_fig(fig, "09_CD20_Aggregates", 'Fig9_Aggregate_Count')

print("\n" + "=" * 60)
print("Wrote 1 figure (CD20+ aggregate count), np.random.seed(42)")
print("  Deterministic count metric — no spatial-distance computation.")
print("  Statistics: Kruskal-Wallis + two-sided Mann-Whitney + Bonferroni.")
print("=" * 60)
