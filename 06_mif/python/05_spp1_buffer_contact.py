#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF step 5 - SPP1+ macrophage buffer neighbourhood in the tumour core
#
# Purpose   : For every SPP1+CD68+ cell, count the unique PD-1+ T cells whose
#             centroid lies within r of it, and compute the contact fraction: the
#             share of macrophages with at least one such T cell within r.
#
# Inputs    : <OUTPUT_DATA>/p6_data.pkl
#
# Outputs   : <OUTPUT_DATA>/Spp1Mac_buffer_contact_TC.csv
#             columns {CD4,CD8,T}_r{20,50,100} plus ..._sampleuniq, ..._frac,
#             ..._per_source and n_spp1_mac; the figures plot the raw union count
#             ..._r{r}, ..._per_source divides it by the number of source cells,
#             ..._frac is the share of source cells with at least one target
#             <OUTPUT_DATA>/Spp1Mac_buffer_Tumor_vendorcheck.csv
#             the tumour-region twin, CSV only, for the vendor cross-check
#             <FIGS_DIR>/10_SPP1Mac_Buffer_Contact/Fig_New_*.{pdf,png}
#
# Usage     : PROJECT_ROOT=/path/to/data python 05_spp1_buffer_contact.py
#             MIF_EXCL=1 PROJECT_ROOT=/path/to/data python 05_spp1_buffer_contact.py
#             The second form writes the exclude-self variant tables.
# -----------------------------------------------------------------------------

import os, pickle, warnings
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import figstyle as fs
warnings.filterwarnings('ignore')

np.random.seed(42)  # reproducibility (scatter jitter)

PX_TO_UM = 0.3455006
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
# MIF_EXCL=1 excludes double-labelled source rows from their own query.
EXCLUDE_SELF = os.environ.get('MIF_EXCL', '0') == '1'
if EXCLUDE_SELF:
    print("MODE: EXCLUDE_SELF — double-labelled source rows excluded from own query", flush=True)
FIGS_DIR = os.path.join(PROJECT_ROOT, "mIF_Validation_Results",
                        "figures" + ("_excl" if EXCLUDE_SELF else ""))
SECTION = "10_SPP1Mac_Buffer_Contact"
SPP1_COL = 'SPP1+CD68+'
TARGETS = {
    'CD4': ('PD-1+CD4+ T', 'PD1+CD4+'),
    'CD8': ('PD-1+CD8+ T', 'PD1+CD8+'),
    'T':   ('PD-1+ T', None),  # union of PD1+CD4+ / PD1+CD8+, double positives once
}
RADII = [20, 50, 100]          # 100 um is the primary radius; 20/50 give the distance series
PRIMARY_R = 100


def get_points_mask_um(df, mask):
    """Centroids (um) of cells where `mask` is True. Returns (0,2) array if none."""
    m = mask
    if not m.any():
        return np.empty((0, 2), dtype=np.float32)
    x = ((df.loc[m, 'XMin'].values + df.loc[m, 'XMax'].values) / 2).astype(np.float32) * PX_TO_UM
    y = ((df.loc[m, 'YMin'].values + df.loc[m, 'YMax'].values) / 2).astype(np.float32) * PX_TO_UM
    return np.column_stack([x, y])


def target_mask(df, tname):
    col = TARGETS[tname][1]
    if col is not None:
        return (df[col] == 1).values
    return ((df['PD1+CD4+'] == 1).values | (df['PD1+CD8+'] == 1).values)


def buffer_total(df, src_mask, tgt_mask, r_um):
    """Number of DISTINCT target cells within r um of ANY source cell
    (union over sources — each target cell counted at most once)."""
    src = get_points_mask_um(df, src_mask)
    tgt = get_points_mask_um(df, tgt_mask)
    if len(src) == 0 or len(tgt) == 0:
        return 0.0
    tree = KDTree(tgt)
    if EXCLUDE_SELF:
        # Per-source exact self-exclusion: a double-labelled source's own
        # target entry is dropped from ITS ball only (another source within
        # radius still keeps it in the union).
        src_dfidx = np.where(src_mask)[0]
        tgt_dfidx = np.where(tgt_mask)[0]
        tgt_pos = {int(dfi): j for j, dfi in enumerate(tgt_dfidx)}
        double = np.logical_and(src_mask, tgt_mask)[src_dfidx]  # aligned to src order
        hit = set()
        for i, p in enumerate(src):
            ball = tree.query_ball_point(p, r_um)
            if double[i]:
                self_j = tgt_pos[int(src_dfidx[i])]
                if self_j in ball:
                    ball.remove(self_j)
            hit.update(ball)
        return float(len(hit))
    hit = set()
    for p in src:
        hit.update(tree.query_ball_point(p, r_um))
    return float(len(hit))


def contact_frac(df, src_mask, tgt_mask, r_um):
    """Fraction of sources with >=1 target cell within r um. NaN if no sources."""
    src = get_points_mask_um(df, src_mask)
    tgt = get_points_mask_um(df, tgt_mask)
    if len(src) == 0:
        return np.nan
    if len(tgt) == 0:
        return 0.0
    tree = KDTree(tgt)
    if EXCLUDE_SELF:
        double_flags = np.logical_and(src_mask, tgt_mask)[src_mask]  # aligned to src order
        hits = 0
        for i, p in enumerate(src):
            nball = len(tree.query_ball_point(p, r_um))
            if double_flags[i]:
                nball -= 1  # self is always inside its own ball
            if nball >= 1:
                hits += 1
        return hits / len(src)
    hits = sum(1 for p in src if len(tree.query_ball_point(p, r_um)) >= 1)
    return hits / len(src)


def compute_sample(df, src_mask):
    """All buffer/contact metrics for one sample (TC or Tumor)."""
    out = {}
    n_src = int(src_mask.sum())
    for tname in TARGETS:
        tm = target_mask(df, tname)
        for r in RADII:
            buf = buffer_total(df, src_mask, tm, r)  # union — each target counted once
            out[f'{tname}_r{r}'] = buf
            out[f'{tname}_r{r}_sampleuniq'] = buf  # alias of the union metric
            # Union count divided by the number of source cells; the figures
            # plot the raw union count above.
            out[f'{tname}_r{r}_per_source'] = (buf / n_src) if n_src else np.nan
            out[f'{tname}_r{r}_frac'] = contact_frac(df, src_mask, tm, r)
    out['n_spp1_mac'] = int(src_mask.sum())
    return out


# ============================================================
# Plotting helpers (identical style to 18/19)
# ============================================================
def save_fig(fig, figname):
    d = os.path.join(FIGS_DIR, SECTION)
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


def r_boxplot_spatial(ax, df, col, title, ylabel):
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
        return False
    fs.panel(ax, dd, present, style='box', title=title, ylabel=ylabel)
    return True


# ============================================================
# Load + compute
# ============================================================
print("Loading p6 data ...", flush=True)
with open(os.path.join(DATA_DIR, 'p6_data.pkl'), 'rb') as f:
    p6_data = pickle.load(f)
print("Loaded.", flush=True)

print(f"\nComputing SPP1+CD68+ buffer/contact (TC, radii {RADII} um) ...", flush=True)
tc_rows = []
for sid, sdata in p6_data['TC'].items():
    df = sdata['df']
    src_mask = (df[SPP1_COL] == 1).values
    row = {'sample_id': sid, 'group': sdata['group']}
    row.update(compute_sample(df, src_mask))
    tc_rows.append(row)
tc_buf = pd.DataFrame(tc_rows)
tc_buf.to_csv(os.path.join(DATA_DIR, 'Spp1Mac_buffer_contact_TC' + ('_excl' if EXCLUDE_SELF else '') + '.csv'), index=False)
print(f"  TC: {len(tc_buf)} samples", flush=True)

print("Computing Tumor-region twin (CSV only, for vendor proximity cross-check) ...", flush=True)
tumor_rows = []
for sid, sdata in p6_data['Tumor'].items():
    df = sdata['df']
    src_mask = (df[SPP1_COL] == 1).values
    row = {'sample_id': sid, 'group': sdata['group']}
    row.update(compute_sample(df, src_mask))
    tumor_rows.append(row)
tumor_buf = pd.DataFrame(tumor_rows)
tumor_buf.to_csv(os.path.join(DATA_DIR, 'Spp1Mac_buffer_Tumor_vendorcheck' + ('_excl' if EXCLUDE_SELF else '') + '.csv'), index=False)
print(f"  Tumor: {len(tumor_buf)} samples", flush=True)

# ============================================================
# Figures
# ============================================================
n_fig = 0

# A. r100 total counts (the primary metric)
print("\nA. 100-um total-count boxplots ...", flush=True)
for tname, (tlabel, _) in TARGETS.items():
    col = f'{tname}_r{PRIMARY_R}'
    fig, ax = plt.subplots(figsize=fs.FIG_PANEL, facecolor='white')
    if r_boxplot_spatial(ax, tc_buf, col,
                         f'{tlabel} within {PRIMARY_R} um of SPP1+ Macrophages (TC)',
                         'Total Cell Count'):
        save_fig(fig, f'Fig_New_{tname}_buffer_r{PRIMARY_R}')
        n_fig += 1
    else:
        plt.close(fig)

# B. contact fraction at r100
print("B. Contact-fraction boxplots ...", flush=True)
for tname, (tlabel, _) in TARGETS.items():
    col = f'{tname}_r{PRIMARY_R}_frac'
    if tc_buf[col].notna().sum() > 6:
        fig, ax = plt.subplots(figsize=fs.FIG_PANEL, facecolor='white')
        if r_boxplot_spatial(ax, tc_buf, col,
                             f'SPP1+ Mac with >=1 {tlabel} within {PRIMARY_R} um (TC)',
                             'Contact Fraction'):
            save_fig(fig, f'Fig_New_{tname}_contact_r{PRIMARY_R}')
            n_fig += 1
        else:
            plt.close(fig)

# C. distance series (r = 20/50/100, mean +/- SEM per group)
print("C. Distance series ...", flush=True)
fig, axes = plt.subplots(1, 3, figsize=(15, 5), facecolor='white')
for ax, (tname, (tlabel, _)) in zip(axes, TARGETS.items()):
    for g in GROUP_ORDER:
        sub = tc_buf[tc_buf['group'] == g]
        if len(sub) == 0:
            continue
        means = [sub[f'{tname}_r{r}'].mean() for r in RADII]
        sems = [sub[f'{tname}_r{r}'].std() / np.sqrt(len(sub)) for r in RADII]
        ax.plot(RADII, means, 'o-', color=USER_COLORS[g], lw=1.5, label=GROUP_LABELS[g])
        ax.fill_between(RADII, [m - s for m, s in zip(means, sems)],
                        [m + s for m, s in zip(means, sems)],
                        color=USER_COLORS[g], alpha=0.15)
    ax.set_xlabel('Radius (um)', fontweight='bold')
    ax.set_ylabel('Mean Total Cell Count', fontweight='bold')
    ax.set_title(f'{tlabel} around SPP1+ Mac (TC)', fontweight='bold')
    ax.set_xticks(RADII)
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.3); ax.set_axisbelow(True)
fig.tight_layout()
save_fig(fig, 'Fig_New_Distance_Series')
n_fig += 1

# ============================================================
# Console summary (per-group medians of the primary metrics)
# ============================================================
print("\n" + "=" * 66)
print("Per-group medians (r100 totals, TC):")
for tname, (tlabel, _) in TARGETS.items():
    col = f'{tname}_r{PRIMARY_R}'
    meds = {g: tc_buf[tc_buf['group'] == g][col].median() for g in GROUP_ORDER}
    print(f"  {tlabel:16s} " + "  ".join(f"{g}: {meds[g]:.1f}" for g in GROUP_ORDER))
print(f"Samples with zero SPP1+CD68+ cells in TC: "
      f"{int((tc_buf['n_spp1_mac'] == 0).sum())}/{len(tc_buf)}")
print(f"Figures saved: {n_fig} -> {os.path.join(FIGS_DIR, SECTION)}")
print("=" * 66)
