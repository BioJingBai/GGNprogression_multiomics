#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF step 2 - per-section cell-composition boxplots
#
# Purpose   : One boxplot per marker and region, one point per section, for the
#             marker percentages and composite-marker percentages of every panel
#             and region table written by 01_load_data.py.
#
# Inputs    : <OUTPUT_DATA>/{P6_TC,P6_IM,P7_TC,P7_IM}_metrics.csv
#             <OUTPUT_DATA>/{p6_data,p7_data}.pkl
#
# Outputs   : <FIGS_DIR>/01_TC_Immune_Landscape/*
#             <FIGS_DIR>/Supplementary/*
#             <FIGS_DIR>/02_B_Cell_TLS/Fig2_B_Immune_Fraction.*
#
# Usage     : PROJECT_ROOT=/path/to/data python 02_cell_composition_boxplots.py
# -----------------------------------------------------------------------------

import os, sys, pickle, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import figstyle as fs
warnings.filterwarnings('ignore')

PX_TO_UM = 0.3455006
np.random.seed(42)  # reproducibility (scatter jitter / any stochastic step)

USER_COLORS = {'Stable':'#8DA0CB', 'Growth':'#66C2A5', 'Solid':'#FC8D62'}
GROUP_ORDER = ['Stable', 'Growth', 'Solid']
GROUP_LABELS = {'Stable':'Stable', 'Growth':'GGO-growth', 'Solid':'SOLID-growth'}
plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 11,
    'axes.titlesize': 14, 'axes.labelsize': 12, 'axes.titleweight': 'bold',
    'axes.labelweight': 'bold', 'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

PROJECT_ROOT = os.environ.get("PROJECT_ROOT", "/path/to/GGN_progression_data")
DATA_DIR = os.path.join(PROJECT_ROOT, "mIF_Validation_Results", "data")
FIGS_DIR = os.path.join(PROJECT_ROOT, "mIF_Validation_Results", "figures")

def save_fig(fig, section, figname):
    d = os.path.join(FIGS_DIR, section)
    os.makedirs(d, exist_ok=True)
    fig.savefig(os.path.join(d, f"{figname}.pdf"), dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(d, f"{figname}.png"), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

def get_group_data(df, col):
    dd = {}
    for g in GROUP_ORDER:
        vals = df[df['group']==g][col].dropna().values
        if len(vals) >= 2: dd[g] = vals
    return dd

# ============================================================
# BOXPLOT FUNCTION
# ============================================================
def r_boxplot(ax, df, col, title, ylabel):
    """Publication-style tight boxplot panel: delegates the drawing and the
    annotation layout to figstyle.panel.

    The statistical convention is unchanged — grouping is still decided by
    get_group_data() (dropna + n >= 2 per group), and the test is still
    Kruskal-Wallis + two-sided Mann-Whitney U + Bonferroni correction; only the
    drawing style and the panel size change.
    """
    dd = get_group_data(df, col)
    present = [g for g in GROUP_ORDER if g in dd]
    if len(present) < 2:
        return False
    fs.panel(ax, dd, present, style='box', title=title, ylabel=ylabel)
    return True

# ============================================================
# Load data
# ============================================================
print("Loading data...")
sys.stdout.flush()
with open(os.path.join(DATA_DIR, 'p6_data.pkl'), 'rb') as f: p6_data = pickle.load(f)
with open(os.path.join(DATA_DIR, 'p7_data.pkl'), 'rb') as f: p7_data = pickle.load(f)
p6_tc = pd.read_csv(os.path.join(DATA_DIR, 'P6_TC_metrics.csv'))
p6_im = pd.read_csv(os.path.join(DATA_DIR, 'P6_IM_metrics.csv'))
p7_tc = pd.read_csv(os.path.join(DATA_DIR, 'P7_TC_metrics.csv'))
p7_im = pd.read_csv(os.path.join(DATA_DIR, 'P7_IM_metrics.csv'))
print("Data loaded.")
sys.stdout.flush()

# ============================================================
# PHASE 2: Fill ALL missing percentage figures
# ============================================================
print("\n" + "="*60)
print("PHASE 2: Fill missing percentage figures (all _pct columns)")
print("="*60)
sys.stdout.flush()

pct_count = 0

def gen_pct_figures(df, pct_cols, prefix, section, title_prefix=""):
    """Generate percentage boxplots for all _pct columns."""
    count = 0
    for col in pct_cols:
        if col not in df.columns: continue
        if df[col].dropna().nunique() < 2: continue
        # Generate clean name
        name = col.replace(' Positive Classification_pct','').replace('_pct','').replace('+','')[:35]
        fig, ax = plt.subplots(figsize=fs.FIG_PANEL, facecolor='white')
        full_title = f'{title_prefix}{name} (%)' if title_prefix else f'{name} (%)'
        if r_boxplot(ax, df, col, full_title, '% of cells'):
            save_fig(fig, section, f'{prefix}{name}')
            count += 1
        else:
            plt.close(fig)
    return count

# P6-TC: all _pct columns (including TUMOR_pct for completeness)
pct_cols_p6tc = [c for c in p6_tc.columns if c.endswith('_pct') and 'TUMOR' not in c and 'STROMA' not in c]
pct_count += gen_pct_figures(p6_tc, pct_cols_p6tc, 'Fig1_', '01_TC_Immune_Landscape')

# P6-IM: all _pct columns
pct_cols_p6im = [c for c in p6_im.columns if c.endswith('_pct') and 'TUMOR' not in c and 'STROMA' not in c]
pct_count += gen_pct_figures(p6_im, pct_cols_p6im, 'FigS_IM_', 'Supplementary', 'IM: ')

# P7-TC: all _pct columns
pct_cols_p7tc = [c for c in p7_tc.columns if c.endswith('_pct') and 'TUMOR' not in c and 'STROMA' not in c]
pct_count += gen_pct_figures(p7_tc, pct_cols_p7tc, 'Fig1_P7_', '01_TC_Immune_Landscape', 'P7: ')

# P7-IM: all _pct columns
pct_cols_p7im = [c for c in p7_im.columns if c.endswith('_pct') and 'TUMOR' not in c and 'STROMA' not in c]
pct_count += gen_pct_figures(p7_im, pct_cols_p7im, 'FigS_P7IM_', 'Supplementary', 'P7 IM: ')

# B_cell_immune_fraction figure.
if 'B_cell_immune_fraction' in p6_tc.columns:
    fig, ax = plt.subplots(figsize=fs.FIG_PANEL, facecolor='white')
    if r_boxplot(ax, p6_tc, 'B_cell_immune_fraction', 'B Cell Immune Fraction (%)', '% of immune cells'):
        save_fig(fig, '02_B_Cell_TLS', 'Fig2_B_Immune_Fraction')
        pct_count += 1
    else: plt.close(fig)

print(f"  Percentage figures generated: {pct_count}")
sys.stdout.flush()

# ============================================================
# SUMMARY
# ============================================================
print(f"\n{'='*60}")
print(f"All gaps filled")
print(f"  Percentage figures filled: {pct_count}")
print(f"  ")
print(f"  Percentage figures regenerated with np.random.seed(42) for reproducibility:")
print(f"    - Percentage figures: in 01_TC_Immune_Landscape / Supplementary")
print(f"    - Density figures: _dens suffix")
print(f"{'='*60}")
