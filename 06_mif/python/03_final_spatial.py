#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF step 3 - nearest-neighbour distances and stromal buffer counts
#
# Purpose   : For every marker pair, compute the per-section distance from each
#             source cell to the nearest target cell (mean and median), and the
#             aSMA+ stromal buffer-zone counts: the number of distinct immune
#             cells within r of any aSMA+ cell.
#
# Inputs    : <OUTPUT_DATA>/p6_data.pkl, <OUTPUT_DATA>/p7_data.pkl
#
# Outputs   : <OUTPUT_DATA>/spatial_results/<pair>_NND.csv          (mean)
#             <OUTPUT_DATA>/spatial_results/<pair>_NND_median.csv   (median)
#             <OUTPUT_DATA>/aSMA_buffer_IM_core4.csv
#             <FIGS_DIR>/<section>/<figure>.{pdf,png}               (NND and buffer panels)
#
# Usage     : PROJECT_ROOT=/path/to/data python 03_final_spatial.py
#             MIF_EXCL=1 PROJECT_ROOT=/path/to/data python 03_final_spatial.py
#             The second form writes the exclude-self variant tables to
#             data/spatial_results_excl/ and figures_excl/.
# -----------------------------------------------------------------------------

import os, pickle, warnings, gc
import numpy as np
import pandas as pd
from scipy.stats import kruskal
from scipy.spatial import KDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import figstyle as fs
warnings.filterwarnings('ignore')

PX_TO_UM = 0.3455006
np.random.seed(42)  # reproducibility

USER_COLORS = {'Stable': '#8DA0CB', 'Growth': '#66C2A5', 'Solid': '#FC8D62'}
GROUP_ORDER = ['Stable', 'Growth', 'Solid']
GROUP_LABELS = {'Stable': 'Stable', 'Growth': 'GGO-growth', 'Solid': 'SOLID-growth'}
plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 11,
    'axes.titlesize': 14, 'axes.labelsize': 12, 'axes.titleweight': 'bold',
    'axes.labelweight': 'bold', 'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

PROJECT_ROOT = os.environ.get("PROJECT_ROOT", "/path/to/GGN_progression_data")
DATA_DIR = os.path.join(PROJECT_ROOT, "mIF_Validation_Results", "data")
# MIF_EXCL=1 excludes double-labelled rows (vendor artifact, e.g. CD20+&CD3+ in
# P7) from their own nearest-neighbour and buffer queries; raw mode is the default.
EXCLUDE_SELF = os.environ.get('MIF_EXCL', '0') == '1'
if EXCLUDE_SELF:
    print("MODE: EXCLUDE_SELF — double-labelled source rows excluded from own query", flush=True)
SPATIAL_DIR = os.path.join(DATA_DIR, "spatial_results" + ("_excl" if EXCLUDE_SELF else ""))
FIGS_DIR = os.path.join(PROJECT_ROOT, "mIF_Validation_Results", "figures" + ("_excl" if EXCLUDE_SELF else ""))
os.makedirs(SPATIAL_DIR, exist_ok=True)

# ============================================================
# Spatial primitives
# ============================================================
def get_points_um(df, marker_col, max_points=None):
    """Centroids (μm) of cells positive for `marker_col`.

    max_points=None keeps ALL cells; a positive value subsamples at random
    (used only for the SOURCE population, where a random sample is an unbiased
    estimator of the population mean/median distance)."""
    mask = (df[marker_col] == 1).values
    n = int(mask.sum())
    if n < 2:
        return None
    if max_points is not None and n > max_points:
        idx = np.random.choice(np.where(mask)[0], max_points, replace=False)
        m = np.zeros(len(df), dtype=bool); m[idx] = True
    else:
        m = mask
    x = ((df.loc[m, 'XMin'].values + df.loc[m, 'XMax'].values) / 2).astype(np.float32) * PX_TO_UM
    y = ((df.loc[m, 'YMin'].values + df.loc[m, 'YMax'].values) / 2).astype(np.float32) * PX_TO_UM
    return np.column_stack([x, y])

def nnd(df, from_col, to_col, agg='mean', max_from=5000):
    """Nearest-neighbour distance (μm) from each `from_col` cell to the nearest
    `to_col` cell. Source subsampled to max_from; target kept in full.

    EXCLUDE_SELF mode: a source row that is also `to_col`-positive is a
    double-labelled vendor-artifact row; its own coordinates are excluded
    from its query (k=2 query — the self hit is the k=1 neighbour)."""
    src_mask = (df[from_col] == 1).values
    tgt_mask = (df[to_col] == 1).values
    src_idx = np.where(src_mask)[0]
    tgt_idx = np.where(tgt_mask)[0]
    if len(src_idx) < 2 or len(tgt_idx) == 0:
        return np.nan
    if max_from is not None and len(src_idx) > max_from:
        src_idx = np.random.choice(src_idx, max_from, replace=False)
    sx = ((df['XMin'].values[src_idx] + df['XMax'].values[src_idx]) / 2).astype(np.float32) * PX_TO_UM
    sy = ((df['YMin'].values[src_idx] + df['YMax'].values[src_idx]) / 2).astype(np.float32) * PX_TO_UM
    src_pts = np.column_stack([sx, sy])
    tx = ((df['XMin'].values[tgt_idx] + df['XMax'].values[tgt_idx]) / 2).astype(np.float32) * PX_TO_UM
    ty = ((df['YMin'].values[tgt_idx] + df['YMax'].values[tgt_idx]) / 2).astype(np.float32) * PX_TO_UM
    tgt_pts = np.column_stack([tx, ty])
    tree = KDTree(tgt_pts)
    if EXCLUDE_SELF:
        double = src_mask[src_idx] & tgt_mask[src_idx]
        if len(tgt_idx) == 1:
            # single target: a double-labelled source IS that target -> no valid neighbour
            dist, _ = tree.query(src_pts, k=1)
            d = np.where(double, np.nan, dist)
        else:
            dist, _ = tree.query(src_pts, k=2)
            d = np.where(double, dist[:, 1], dist[:, 0])
        d = d[np.isfinite(d)]
        if len(d) == 0:
            return np.nan
        return {'mean': np.mean, 'median': np.median, 'min': np.min}[agg](d)
    dist, _ = tree.query(src_pts, k=1)
    return {'mean': np.mean, 'median': np.median, 'min': np.min}[agg](dist)

def compute_buffer(df, asma_col, immune_col, radius_um):
    """Number of DISTINCT immune cells within `radius_um` of ANY aSMA cell
    (union over sources — each target cell counted at most once). Both
    populations kept in full.

    EXCLUDE_SELF mode: a double-labelled source row (positive for both
    markers) hits ITSELF at distance 0 — that self-hit is removed from its
    own ball only (another source within radius still keeps the cell)."""
    asma_pts = get_points_um(df, asma_col, max_points=None)
    tgt_mask = (df[immune_col] == 1).values
    tgt_pos = np.where(tgt_mask)[0]
    if asma_pts is None or len(tgt_pos) < 2:
        return np.nan
    x = ((df.iloc[tgt_pos]['XMin'].values + df.iloc[tgt_pos]['XMax'].values) / 2).astype(np.float32) * PX_TO_UM
    y = ((df.iloc[tgt_pos]['YMin'].values + df.iloc[tgt_pos]['YMax'].values) / 2).astype(np.float32) * PX_TO_UM
    tree = KDTree(np.column_stack([x, y]))
    hit = set()
    if EXCLUDE_SELF:
        dbl = ((df[asma_col] == 1) & (df[immune_col] == 1)).values
        src_dfidx = np.where((df[asma_col] == 1).values)[0]
        tgt_pos_map = {int(dfi): j for j, dfi in enumerate(tgt_pos)}
        for k, pt in enumerate(asma_pts):
            ball = tree.query_ball_point(pt, radius_um)
            if dbl[src_dfidx[k]]:
                self_j = tgt_pos_map[int(src_dfidx[k])]
                if self_j in ball:
                    ball.remove(self_j)
            hit.update(ball)
    else:
        for pt in asma_pts:
            hit.update(tree.query_ball_point(pt, radius_um))
    return float(len(hit))

# ============================================================
# Load data
# ============================================================
print("Loading data (large pkl files)...", flush=True)
with open(os.path.join(DATA_DIR, 'p6_data.pkl'), 'rb') as f:
    p6_data = pickle.load(f)
with open(os.path.join(DATA_DIR, 'p7_data.pkl'), 'rb') as f:
    p7_data = pickle.load(f)
print("Data loaded.", flush=True)

# ============================================================
# PHASE 1 — recompute every NND (mean + median)
# ============================================================
CONFIGS = [
    # B-cell interactions (P6-TC)
    ('p6_TC', p6_data['TC'], 'CD20 Positive Classification', 'CD8 Positive Classification', 'B_CD8_NND'),
    ('p6_TC', p6_data['TC'], 'CD20 Positive Classification', 'CD4 Positive Classification', 'B_CD4_NND'),
    ('p6_TC', p6_data['TC'], 'CD20 Positive Classification', 'GZMB+CD8+', 'B_GZMB_CD8_NND'),
    ('p6_TC', p6_data['TC'], 'CD20 Positive Classification', 'PD1+CD8+', 'B_PD1_CD8_NND'),
    ('p6_TC', p6_data['TC'], 'CD20 Positive Classification', 'PD1+CD4+', 'B_PD1_CD4_NND'),
    ('p6_TC', p6_data['TC'], 'CD20 Positive Classification', 'PD1 Positive Classification', 'B_PD1_NND'),
    ('p6_TC', p6_data['TC'], 'CD20 Positive Classification', 'CD68 Positive Classification', 'B_CD68_NND'),
    # SPP1+ macrophage interactions (P6-TC)
    ('p6_TC', p6_data['TC'], 'SPP1+CD68+', 'CD8 Positive Classification', 'SPP1Mac_CD8_NND'),
    ('p6_TC', p6_data['TC'], 'SPP1+CD68+', 'PD1+CD8+', 'SPP1Mac_PD1_CD8_NND'),
    ('p6_TC', p6_data['TC'], 'SPP1+CD68+', 'CD4 Positive Classification', 'SPP1Mac_CD4_NND'),
    ('p6_TC', p6_data['TC'], 'SPP1+CD68+', 'PD1+CD4+', 'SPP1Mac_PD1_CD4_NND'),
    ('p6_TC', p6_data['TC'], 'SPP1+CD68+', 'PD1 Positive Classification', 'SPP1Mac_PD1_NND'),
    ('p6_TC', p6_data['TC'], 'SPP1+CD68+', 'CD20 Positive Classification', 'SPP1Mac_CD20_NND'),
    ('p6_TC', p6_data['TC'], 'SPP1+CD68+', 'aSMA Positive Classification', 'SPP1Mac_aSMA_NND'),
    # aSMA barrier (P6-IM and P6-TC)
    ('p6_IM', p6_data['IM'], 'aSMA Positive Classification', 'CD8 Positive Classification', 'aSMA_CD8_IM_NND'),
    ('p6_IM', p6_data['IM'], 'aSMA Positive Classification', 'CD20 Positive Classification', 'aSMA_CD20_IM_NND'),
    ('p6_IM', p6_data['IM'], 'aSMA Positive Classification', 'CD4 Positive Classification', 'aSMA_CD4_IM_NND'),
    ('p6_IM', p6_data['IM'], 'aSMA Positive Classification', 'SPP1+CD68+', 'aSMA_SPP1Mac_IM_NND'),
    ('p6_TC', p6_data['TC'], 'aSMA Positive Classification', 'CD8 Positive Classification', 'aSMA_CD8_TC_NND'),
    ('p6_TC', p6_data['TC'], 'aSMA Positive Classification', 'CD20 Positive Classification', 'aSMA_CD20_TC_NND'),
    # HLA-DR+ macrophage interactions (P7)
    ('p7_TC', p7_data['TC'], 'HLA-DR+CD68+', 'CD3 Positive Classification', 'HLADRMac_CD3_NND'),
    ('p7_TC', p7_data['TC'], 'HLA-DR+CD68+', 'CD20 Positive Classification', 'HLADRMac_CD20_NND'),
    ('p7_IM', p7_data['IM'], 'HLA-DR+CD68+', 'PANCK Positive Classification', 'HLADRMac_PANCK_IM_NND'),
    ('p7_TC', p7_data['TC'], 'HLA-DR+CD68+', 'PANCK Positive Classification', 'HLADRMac_PANCK_TC_NND'),
    # CCR7-CCL19 (P7-TC)
    ('p7_TC', p7_data['TC'], 'CCR7 Positive Classification', 'CCL19 Positive Classification', 'CCR7_CCL19_NND'),
]

print("\n" + "=" * 60)
print("PHASE 1: nearest-neighbour distances (mean / median / min)")
print("=" * 60)
for panel, data_dict, from_col, to_col, fname in CONFIGS:
    # Mean, median and minimum are each written out; the panels display the median.
    for agg, suffix in [('mean', ''), ('median', '_median'), ('min', '_min')]:
        rows = []
        for sid, sdata in data_dict.items():
            v = nnd(sdata['df'], from_col, to_col, agg=agg)
            if not np.isnan(v):
                rows.append({'sample_id': sid, 'group': sdata['group'], 'obs_nnd': v})
        df = pd.DataFrame(rows)
        if len(df) > 5:
            df.to_csv(os.path.join(SPATIAL_DIR, f"{fname}{suffix}.csv"), index=False)
    gc.collect()
    print(f"  {fname}: n={len(df)}", flush=True)

# B-NND (B-cell self-clustering; ALL B cells, k=2 excludes self)
print("  B-NND (self) ...", flush=True)
for agg, suffix in [('mean', ''), ('median', '_median')]:
    rows = []
    for sid, sdata in p6_data['TC'].items():
        pts = get_points_um(sdata['df'], 'CD20 Positive Classification', max_points=None)
        if pts is not None and len(pts) >= 3:
            tree = KDTree(pts)
            dist, _ = tree.query(pts, k=2)
            d = dist[:, 1]  # nearest OTHER cell
            rows.append({'sample_id': sid, 'group': sdata['group'],
                         'obs_nnd': np.mean(d) if agg == 'mean' else np.median(d)})
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(SPATIAL_DIR, f"B_NND{suffix}.csv"), index=False)
        print(f"  B_NND{suffix}: n={len(rows)}", flush=True)

# aSMA -> PANCK (barrier-to-tumour distance; target PANCK kept in full)
for region, data_dict in [('IM', p6_data['IM']), ('TC', p6_data['TC'])]:
    for agg, suffix in [('mean', ''), ('median', '_median')]:
        rows = []
        for sid, sdata in data_dict.items():
            v = nnd(sdata['df'], 'aSMA Positive Classification', 'PANCK Positive Classification', agg=agg)
            if not np.isnan(v):
                rows.append({'sample_id': sid, 'group': sdata['group'], 'obs_nnd': v})
        if rows:
            pd.DataFrame(rows).to_csv(os.path.join(SPATIAL_DIR, f"aSMA_to_PANCK_{region}{suffix}.csv"), index=False)
    print(f"  aSMA_to_PANCK_{region} (mean + median)", flush=True)

# ============================================================
# PHASE 2 — aSMA buffer zone (vendor-matching metric: DISTINCT immune cells
# within r of ANY aSMA+ cell; both populations kept in full)
# ============================================================
print("\n" + "=" * 60)
print("PHASE 2: aSMA buffer zone (correct method)")
print("=" * 60)
RADII = [20, 50]  # μm
IMMUNE_BUF = {
    'CD8': 'CD8 Positive Classification',
    'CD4': 'CD4 Positive Classification',
    'CD20': 'CD20 Positive Classification',
    'CD68': 'CD68 Positive Classification',
}
buffer_rows = []
for sid, sdata in p6_data['IM'].items():
    df = sdata['df']
    row = {'sample_id': sid, 'group': sdata['group']}
    for name, col in IMMUNE_BUF.items():
        for r in RADII:
            row[f'{name}_r{r}'] = compute_buffer(df, 'aSMA Positive Classification', col, r)
    buffer_rows.append(row)
buffer_df = pd.DataFrame(buffer_rows)
# Core-4 buffer table (CD8/CD4/CD20/CD68 x r20/r50).
buffer_df.to_csv(os.path.join(DATA_DIR, 'aSMA_buffer_IM_core4' + ('_excl' if EXCLUDE_SELF else '') + '.csv'), index=False)
print(f"  Buffer data (core 4 markers, r20/r50): {len(buffer_df)} samples", flush=True)

# ============================================================
# PHASE 3 — figures
# ============================================================
def get_group_data(df, col):
    dd = {}
    for g in GROUP_ORDER:
        vals = df[df['group'] == g][col].dropna().values
        if len(vals) >= 2:
            dd[g] = vals
    return dd

def r_boxplot_spatial(ax, df, col, title, ylabel):
    """Publication-style tight boxplot panel, delegating to figstyle.panel.

    The statistical convention and the values are unchanged: grouping is still
    decided by get_group_data() (dropna + n >= 2 per group), and the test is
    still Kruskal-Wallis + two-sided Mann-Whitney U + Bonferroni correction;
    only the drawing style and the panel size change.
    """
    dd = get_group_data(df, col)
    present = [g for g in GROUP_ORDER if g in dd]
    if len(present) < 2:
        return False
    fs.panel(ax, dd, present, style='box', title=title, ylabel=ylabel)
    return True

def save_fig(fig, section, figname):
    d = os.path.join(FIGS_DIR, section)
    os.makedirs(d, exist_ok=True)
    fig.savefig(os.path.join(d, f"{figname}.pdf"), dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(d, f"{figname}.png"), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)

# (csv base name, title, region, figure base name)
NND_FIG_DEFS = {
    '02_B_Cell_TLS': [
        ('B_NND', 'B Cell Clustering: B→B NND', 'TC', 'Fig2_B_NND'),
        ('B_CD8_NND', 'B→CD8+ T', 'TC', 'Fig2_B_CD8_NND'),
        ('B_CD4_NND', 'B→CD4+ T', 'TC', 'Fig2_B_CD4_NND'),
        ('B_GZMB_CD8_NND', 'B→GZMB+CD8+ T', 'TC', 'Fig2_B_GZMB_CD8_NND'),
        ('B_PD1_CD8_NND', 'B→PD-1+CD8+ T', 'TC', 'Fig2_B_PD1_CD8_NND'),
        ('B_PD1_CD4_NND', 'B→PD-1+CD4+ T', 'TC', 'Fig2_B_PD1_CD4_NND'),
        ('B_PD1_NND', 'B→PD-1+ Cells', 'TC', 'Fig2_B_PD1_NND'),
        ('B_CD68_NND', 'B→CD68+ Mac', 'TC', 'Fig2_B_CD68_NND'),
    ],
    '03_Macrophage_Polarization': [
        ('SPP1Mac_CD8_NND', 'SPP1+Mac→CD8+ T', 'TC', 'Fig3_SPP1Mac_CD8_NND'),
        ('SPP1Mac_PD1_CD8_NND', 'SPP1+Mac→PD-1+CD8+ T', 'TC', 'Fig3_SPP1Mac_PD1_CD8_NND'),
        ('SPP1Mac_CD4_NND', 'SPP1+Mac→CD4+ T', 'TC', 'Fig3_SPP1Mac_CD4_NND'),
        ('SPP1Mac_PD1_CD4_NND', 'SPP1+Mac→PD-1+CD4+ T', 'TC', 'Fig3_SPP1Mac_PD1_CD4_NND'),
        ('SPP1Mac_PD1_NND', 'SPP1+Mac→PD-1+ Cells', 'TC', 'Fig3_SPP1Mac_PD1_NND'),
        ('SPP1Mac_CD20_NND', 'SPP1+Mac→CD20+ B', 'TC', 'Fig3_SPP1Mac_CD20_NND'),
        ('SPP1Mac_aSMA_NND', 'SPP1+Mac→aSMA+ Fibroblast', 'TC', 'Fig3_SPP1Mac_aSMA_NND'),
        ('HLADRMac_CD3_NND', 'HLA-DR+Mac→CD3+ T', 'TC', 'Fig3_HLADRMac_CD3_NND'),
        ('HLADRMac_CD20_NND', 'HLA-DR+Mac→CD20+ B', 'TC', 'Fig3_HLADRMac_CD20_NND'),
        ('HLADRMac_PANCK_IM_NND', 'HLA-DR+Mac→PANCK+ (IM)', 'IM', 'Fig3_HLADRMac_PANCK_IM_NND'),
        ('HLADRMac_PANCK_TC_NND', 'HLA-DR+Mac→PANCK+ (TC)', 'TC', 'Fig3_HLADRMac_PANCK_TC_NND'),
    ],
    '05_aSMA_Barrier': [
        ('aSMA_CD8_IM_NND', 'aSMA+→CD8+ T', 'IM', 'Fig5_aSMA_CD8_IM_NND'),
        ('aSMA_CD20_IM_NND', 'aSMA+→CD20+ B', 'IM', 'Fig5_aSMA_CD20_IM_NND'),
        ('aSMA_CD4_IM_NND', 'aSMA+→CD4+ T', 'IM', 'Fig5_aSMA_CD4_IM_NND'),
        ('aSMA_SPP1Mac_IM_NND', 'aSMA+→SPP1+Mac', 'IM', 'Fig5_aSMA_SPP1Mac_IM_NND'),
        ('aSMA_CD8_TC_NND', 'aSMA+→CD8+ T (TC Ctrl)', 'TC', 'Fig5_aSMA_CD8_TC_NND'),
        ('aSMA_CD20_TC_NND', 'aSMA+→CD20+ B (TC Ctrl)', 'TC', 'Fig5_aSMA_CD20_TC_NND'),
        ('aSMA_to_PANCK_IM', 'aSMA→PANCK Distance (IM)', '', 'Fig5_aSMA_to_PANCK_nnd_IM'),
        ('aSMA_to_PANCK_TC', 'aSMA→PANCK Distance (TC)', '', 'Fig5_aSMA_to_PANCK_nnd_TC'),
    ],
    '07_TLS_Chemokines': [
        ('CCR7_CCL19_NND', 'CCR7+→CCL19+', 'TC', 'Fig7_CCR7_CCL19_NND'),
    ],
}

fig_count = 0
print("\n" + "=" * 60)
print("PHASE 3: generating NND figures (mean + median)")
print("=" * 60)
for section, defs in NND_FIG_DEFS.items():
    for csv_base, title, region, figname in defs:
        # mean (original name)
        csv_path = os.path.join(SPATIAL_DIR, f"{csv_base}.csv")
        if os.path.exists(csv_path):
            fig, ax = plt.subplots(figsize=fs.FIG_PANEL, facecolor='white')
            if r_boxplot_spatial(ax, pd.read_csv(csv_path), 'obs_nnd',
                                 (f'{title} ({region})' if region else title), 'Mean NND (μm)'):
                save_fig(fig, section, figname); fig_count += 1
            else:
                plt.close(fig)
        # median (_median suffix)
        csv_path = os.path.join(SPATIAL_DIR, f"{csv_base}_median.csv")
        if os.path.exists(csv_path):
            fig, ax = plt.subplots(figsize=fs.FIG_PANEL, facecolor='white')
            if r_boxplot_spatial(ax, pd.read_csv(csv_path), 'obs_nnd',
                                 (f'{title} ({region})' if region else title), 'Median NND (μm)'):
                save_fig(fig, section, f"{figname}_median"); fig_count += 1
            else:
                plt.close(fig)
print(f"  NND figures: {fig_count}", flush=True)

# ============================================================
# Buffer figures
# ============================================================
print("  Generating buffer figures...", flush=True)
buf_fig = 0
S5 = '05_aSMA_Barrier'
for immune_name in ['CD8', 'CD4', 'CD20', 'CD68']:
    for r in RADII:
        col = f'{immune_name}_r{r}'
        if col in buffer_df.columns and buffer_df[col].notna().sum() > 6:
            fig, ax = plt.subplots(figsize=fs.FIG_PANEL, facecolor='white')
            if r_boxplot_spatial(ax, buffer_df, col,
                                 f'{immune_name}+ within {r}μm of aSMA+ Cells (IM)', 'Total Cell Count'):
                save_fig(fig, S5, f'Fig5_Buffer_{immune_name}_r{r}'); buf_fig += 1
            else:
                plt.close(fig)

# Multi-panel overview (20 μm)
fig, axes = plt.subplots(2, 2, figsize=fs.FIG_2X2, facecolor='white')
for i, immune_name in enumerate(['CD8', 'CD4', 'CD20', 'CD68']):
    ax = axes[i // 2, i % 2]
    col = f'{immune_name}_r20'
    dd = get_group_data(buffer_df, col)
    gs = [g for g in GROUP_ORDER if g in dd]
    if len(gs) >= 2:
        kw_s, kw_p = kruskal(*[dd[g] for g in gs])
        fs.panel(ax, dd, gs, style='box', brackets=False, kw_in_title=False,
                 title=f'{immune_name}+ | 20$\\mu$m | KW p='
                       f'{"<0.001" if kw_p < 0.001 else f"{kw_p:.3f}"}',
                 ylabel='Total Cell Count')
    else:
        ax.set_ylabel('Total Cell Count', fontsize=fs.FS['axes.labelsize'])
        ax.set_yticks([])
        fs.classic(ax)
fig.suptitle('Immune Cells within 20μm of aSMA+ Cells (IM)', fontsize=8.4)
plt.tight_layout()
save_fig(fig, S5, 'Fig5_Buffer_20um_Overview'); buf_fig += 1

# Line chart: total count vs distance
fig, axes = plt.subplots(1, 4, figsize=(22, 5.5), facecolor='white')
for i, immune_name in enumerate(['CD8', 'CD4', 'CD20', 'CD68']):
    ax = axes[i]
    for g in GROUP_ORDER:
        gdf = buffer_df[buffer_df['group'] == g]
        if len(gdf) >= 3:
            means = [gdf[f'{immune_name}_r{r}'].mean() for r in RADII]
            sems = [gdf[f'{immune_name}_r{r}'].std() / np.sqrt(len(gdf)) for r in RADII]
            ax.errorbar(RADII, means, yerr=sems, color=USER_COLORS[g], marker='o', ms=8,
                        lw=2.5, capsize=5, label=GROUP_LABELS[g])
    ax.set_xlabel('Distance from aSMA+ (μm)', fontsize=10)
    ax.set_ylabel(f'{immune_name}+ Total Count', fontsize=10)
    ax.set_title(immune_name, fontweight='bold', fontsize=12)
    ax.legend(fontsize=8); ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)
fig.suptitle('aSMA+ Buffer Zone: Immune Cell Count by Distance (IM)', fontweight='bold', fontsize=16)
plt.tight_layout()
save_fig(fig, S5, 'Fig5_Buffer_Distance_Lines'); buf_fig += 1
print(f"  Buffer figures: {buf_fig}", flush=True)

# ============================================================
print(f"\n{'=' * 60}")
print(f"Final spatial pipeline complete")
print(f"  NND figures (mean + median): {fig_count}")
print(f"  Buffer figures: {buf_fig}")
print(f"  Coordinates: pixel × {PX_TO_UM} = μm")
print(f"  NND target cells: ALL (no subsampling bias)")
print(f"  Seed: 42 (reproducible)")
print(f"  Statistics: Kruskal-Wallis + Mann-Whitney (two-sided) + Bonferroni")
print(f"{'=' * 60}")
