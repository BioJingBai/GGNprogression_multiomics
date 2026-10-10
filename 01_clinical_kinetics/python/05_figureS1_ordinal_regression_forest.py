#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Stratified aggressiveness gradients and the adjusted forest plot
#
# Purpose   : Draw four panels: a two-axis KDE density plane, the aggressiveness
#             gradient within baseline strata (baseline pGGN and baseline PSN),
#             the effect-size forest plot before and after adjustment for baseline
#             3D CTR and volume, and the symmetric endpoint test within last-CT
#             strata.
#
# Inputs    : ./kinetics.csv (built by 01_build_nodule_kinetics_table.py)
#
# Outputs   : ./Fig13_twoaxis_density.png/.pdf, ./Fig14_baseline_stratified.png/.pdf,
#             ./Fig15_ctr_adjusted_forest.png/.pdf, ./Fig16_finalCT_stratified.png/.pdf
#
# Usage     : PROJECT_ROOT=/path/to/data python3 05_figureS1_ordinal_regression_forest.py
# -----------------------------------------------------------------------------
"""
Four supplementary panels, all built from kinetics.csv in this folder:

  Fig13_twoaxis_density     two-axis plane: per-group KDE contours (50% / 95%) +
                            group centroids + fixed reference lines + per-group n
  Fig14_baseline_stratified aggressiveness gradient within baseline strata
                            (baseline pGGN / baseline PSN)
  Fig15_ctr_adjusted_forest effect-size comparison before and after adjusting for
                            baseline 3D CTR and baseline volume
  Fig16_finalCT_stratified  aggressiveness gradient within last-CT strata, the
                            endpoint counterpart of Fig14

Definitions: primary VDT is set by this folder's MODE; mass doubling time = 3D;
pathology 0 benign < 1 AAH < 2 AIS < 3 MIA < 4 IAC.
Baseline nodule-type code (source table 0/1/2) is binarised as 0 = pGGN and
1-2 = PSN (the heterogeneous ground-glass category folded into PSN).
"""
# -----------------------------------------------------------------------------

import os, csv, json, textwrap
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
import statsmodels.api as sm
from statsmodels.miscmodels.ordinal_model import OrderedModel

OUT = os.path.dirname(os.path.abspath(__file__)) + "/"
MODE = {'VDT2d': 'VDT 2D ellipsoid', 'VDT2d_fill3d': 'VDT 2D-first, 3D filling where 2D unavailable',
        'VDT3d': 'VDT 3D segmentation'}

df = pd.read_csv(OUT + "kinetics.csv")
df['group'] = df['group'].replace({'GK253GO-growth': 'GGO-growth'})
df = df[df['group'].isin(['Stable', 'GGO-growth', 'SOLID-growth'])].copy()
df['pn'] = pd.to_numeric(df['patho'], errors='coerce')
df.loc[df['pn'] < 0, 'pn'] = np.nan
df['IAC'] = (df['pn'] == 4).astype(float)
df.loc[df['pn'].isna(), 'IAC'] = np.nan
CUR_MODE = str(df['MODE'].dropna().iloc[0]) if df['MODE'].notna().any() else 'unknown'

GROUPS = ['Stable', 'GGO-growth', 'SOLID-growth']
COLORS = {'Stable': '#4C72B0', 'GGO-growth': '#55A868', 'SOLID-growth': '#C44E52'}
GRADE_LAB = ['Benign', 'AAH', 'AIS', 'MIA', 'IAC']
GRADE_COL = ['#D9D9D9', '#BDD7E7', '#8FBEE0', '#5B9BD5', '#1F4E79']

rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'],
                 'font.size': 8.5, 'axes.linewidth': 0.7, 'xtick.direction': 'out',
                 'ytick.direction': 'out', 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                 'figure.dpi': 300})


# --------------------------------------------------------------- statistics helpers
def jt_test(groups):
    """Jonckheere-Terpstra test for an increasing trend (standard normal approximation, mid-ranks with tie correction)."""
    k = len(groups)
    n = [len(g) for g in groups]
    N = sum(n)
    U = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            a, b = groups[i], groups[j]
            U += (a[:, None] < b[None, :]).sum() + 0.5 * (a[:, None] == b[None, :]).sum()
    EU = (N * N - sum(x * x for x in n)) / 4.0
    VU = (N * N * (2 * N + 3) - sum(x * x * (2 * x + 3) for x in n)) / 72.0
    z = (U - EU) / np.sqrt(VU)
    return z, float(2 * stats.norm.sf(abs(z))), n


def bh(ps):
    ps = np.asarray(ps, float)
    k = len(ps)
    order = np.argsort(ps)
    adj = np.empty(k)
    adj[order] = np.minimum.accumulate((ps[order] * k / np.arange(1, k + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def star(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def pfmt(p):
    return 'P < 1e-300' if p == 0 else ('P = %.2g' % p)


LOG = []


def log(s=''):
    print(s)
    LOG.append(str(s))


def wrap_caption(txt, width=170):
    """Wrap the caption: avoid a single over-long line stretching the bbox_inches='tight' canvas."""
    return textwrap.fill(txt, width=width)


def panel_letter(ax, L, dx=-0.16, dy=1.06):
    ax.text(dx, dy, L, transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')


# ======================================================================
# Fig 13: two-axis plane (whole-cohort KDE background + three-group median/IQR + ordered progression arrows + marginal ridges)
#   x = growth-rate axis: log10(365.25 / VDT)  (doublings/yr; larger = faster; monotonically inverse to VDT)
#   y = densification axis: annualized CTR change (/yr)
# ======================================================================
CLIP_SG, CLIP_SOLID = 0.2, 0.3
VREF = 400.0                                  # VDT = 400 d (fast-growth landmark)
GRX = lambda vdt: np.log10(365.25 / vdt)      # growth-rate axis transform
xv = GRX(VREF)
XL = (-2.60, 1.05)                            # display window: VDT about 32 - 145000 days
YL = (-CLIP_SG - 0.02, CLIP_SOLID + 0.02)
GRIDX, GRIDY = 220, 220

pts = {}
for g in GROUPS:
    d = df[(df['group'] == g)].dropna(subset=['VDT', 'ann_CTR3d'])
    d = d[d['VDT'] > 0]
    if g == 'SOLID-growth':
        d = d[(d['ann_CTR3d'] >= 0) & (d['ann_CTR3d'] <= CLIP_SOLID)]
    else:
        d = d[d['ann_CTR3d'].abs() <= CLIP_SG]
    pts[g] = (GRX(d['VDT'].values), d['ann_CTR3d'].values)

gx = np.linspace(*XL, GRIDX)
gy = np.linspace(*YL, GRIDY)
GX, GY = np.meshgrid(gx, gy)
POS = np.vstack([GX.ravel(), GY.ravel()])

fig = plt.figure(figsize=(6.4, 5.4))
gs = fig.add_gridspec(2, 2, width_ratios=[5, 1.25], height_ratios=[1.25, 5],
                      wspace=0.07, hspace=0.07)
ax_top = fig.add_subplot(gs[0, 0])
ax = fig.add_subplot(gs[1, 0])
ax_right = fig.add_subplot(gs[1, 1])

# --- Main plane: whole-cohort 2-D KDE grey background (continuum) + group median / IQR cross ---
pool_x = np.concatenate([pts[g][0] for g in GROUPS])
pool_y = np.concatenate([pts[g][1] for g in GROUPS])
Zall = stats.gaussian_kde(np.vstack([pool_x, pool_y]))(POS).reshape(GX.shape)
zsa = np.sort(Zall.ravel())[::-1]
cuma = np.cumsum(zsa) / zsa.sum()
LV = sorted(zsa[np.searchsorted(cuma, q)] for q in (0.25, 0.50, 0.75))
ax.contourf(GX, GY, Zall, levels=[LV[0], Zall.max() * 1.01], colors=['#C9C9C9'],
            alpha=0.45, zorder=1)
ax.contour(GX, GY, Zall, levels=LV, colors=['#8C8C8C'], linewidths=[0.6, 0.9, 0.6],
           linestyles=[':', '-', ':'], zorder=1)

MED, IQR = {}, {}
for g in GROUPS:
    xs, ys = pts[g]
    mx, my = float(np.median(xs)), float(np.median(ys))
    qx = np.percentile(xs, [25, 75]); qy = np.percentile(ys, [25, 75])
    MED[g] = (mx, my); IQR[g] = (qx, qy)
    ax.plot([qx[0], qx[1]], [my, my], lw=1.2, color=COLORS[g], alpha=0.9,
            solid_capstyle='butt', zorder=5)
    ax.plot([mx, mx], [qy[0], qy[1]], lw=1.2, color=COLORS[g], alpha=0.9,
            solid_capstyle='butt', zorder=5)
    ax.scatter([mx], [my], s=120, marker='o', c=COLORS[g], edgecolors='black',
               linewidths=0.9, zorder=6)
    log(f'[Fig13] {g}: n(display)={len(xs)}, median VDT={365.25/10**mx:.0f} d, '
        f'median annCTR={my:.4f}/yr, IQR VDT={365.25/10**qx[0]:.0f}-{365.25/10**qx[1]:.0f} d, '
        f'IQR annCTR={qy[0]:.4f}-{qy[1]:.4f}/yr')

LBL_OFF = {'Stable': (0, 26), 'GGO-growth': (-6, -30), 'SOLID-growth': (26, 24)}
# The first segment (Stable->GGO) advances horizontally along the growth-rate axis;
# the second (GGO->SOLID) has only a few percentage points of vertical separation, so
# drawn straight from the centroids it would be hidden by the two dots and the IQR
# bars; the second segment is therefore shifted right and drawn as a side-by-side
# upward arrow.
ax.annotate('', xy=MED['GGO-growth'], xytext=MED['Stable'],
            arrowprops=dict(arrowstyle='-|>', color='#3A3A3A', lw=1.3,
                            shrinkA=9, shrinkB=11, alpha=0.9, zorder=5))
_g, _s = MED['GGO-growth'], MED['SOLID-growth']
_dx = 0.30 * (XL[1] - XL[0]) / 3.65
ax.annotate('', xy=(_s[0] + _dx, _s[1] - 0.028), xytext=(_g[0] + _dx, _g[1] + 0.010),
            arrowprops=dict(arrowstyle='-|>', color='#3A3A3A', lw=1.3, alpha=0.9,
                            shrinkA=0, shrinkB=0), zorder=5)
for g in GROUPS:
    ax.annotate(f'{g}\nn = {len(pts[g][0])}', MED[g], textcoords='offset points',
                xytext=LBL_OFF[g], ha='center', fontsize=7.6, fontweight='bold',
                color=COLORS[g], zorder=7,
                bbox=dict(boxstyle='round,pad=0.18', fc='white', ec='none', alpha=0.85))

ax.axhline(0, ls='--', lw=0.8, c='grey', zorder=4)
ax.axvline(xv, ls=':', lw=1.0, c='grey', zorder=4)
ax.set_xlim(*XL)
ax.set_ylim(*YL)
XT = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
ax.set_xticks(np.log10(XT))
ax.set_xticklabels([f'{t:g}' for t in XT])
ax.set_xlabel('Growth rate  (doublings per year, log scale;  365 / VDT)', fontsize=8.2)
ax.set_ylabel('Annualized CTR change (/yr)  —  densification axis', fontsize=8.2)
ax.text(xv + 0.03, YL[1] - 0.012, f'VDT = {VREF:.0f} d', rotation=90,
        va='top', ha='left', fontsize=6.4, color='grey')
ax.text(XL[1] - 0.03, -0.055, 'no densification (ΔCTR/yr = 0)', ha='right', va='top',
        fontsize=6.4, color='grey')
ax.text(XL[0] + 0.06, YL[1] - 0.012, 'all nodules, pooled density', ha='left', va='top',
        fontsize=6.4, color='#8C8C8C')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# --- 1-D KDE ridges for each axis (offset by group) ---
xgrid = np.linspace(*XL, 300)
ygrad = np.linspace(*YL, 300)
for gi, g in enumerate(GROUPS):
    xs, ys = pts[g]
    kx = stats.gaussian_kde(xs)(xgrid); kx = kx / kx.max() * 0.62
    ax_top.fill_between(xgrid, gi, gi + kx, color=COLORS[g], alpha=0.55, lw=0)
    ax_top.plot(xgrid, gi + kx, color=COLORS[g], lw=0.8)
    ky = stats.gaussian_kde(ys)(ygrad); ky = ky / ky.max() * 0.62
    ax_right.fill_betweenx(ygrad, gi, gi + ky, color=COLORS[g], alpha=0.55, lw=0)
    ax_right.plot(gi + ky, ygrad, color=COLORS[g], lw=0.8)
ax_top.set_xlim(*XL)
ax_top.set_ylim(-0.05, 3.0)
ax_top.set_yticks([0.5, 1.5, 2.5])
ax_top.set_yticklabels(GROUPS, fontsize=6.6)
ax_top.set_xticks(np.log10(XT))
ax_top.tick_params(axis='x', labelbottom=False)
ax_top.axvline(xv, ls=':', lw=1.0, c='grey')
ax_top.set_title('Two-axis kinetics: the three trajectory groups as ordered anchors of one continuum',
                 fontsize=9.0, pad=6)
for s in ['top', 'right', 'bottom']:
    ax_top.spines[s].set_visible(False)
ax_top.spines['left'].set_linewidth(0.5)

ax_right.set_ylim(*YL)
ax_right.set_xlim(-0.05, 3.0)
ax_right.set_xticks([0.5, 1.5, 2.5])
ax_right.set_xticklabels(GROUPS, fontsize=6.6, rotation=90)
ax_right.tick_params(axis='y', labelleft=False)
ax_right.axhline(0, ls='--', lw=0.8, c='grey')
for s in ['top', 'left', 'right']:
    ax_right.spines[s].set_visible(False)
ax_right.spines['bottom'].set_linewidth(0.5)

CAP13 = (f'Main plane: grey shading, pooled 2-D kernel density of all nodules; coloured cross-bars, group median '
         f'with interquartile range on each axis (both axes time-normalized); arrows, the ordered progression. '
         f'Marginals: group 1-D kernel densities. Dashed line, no densification (ΔCTR/yr = 0); dotted line, '
         f'slow/fast growth reference (VDT = 400 d). Display: |ΔCTR/yr| ≤ {CLIP_SG} (Stable, GGO-growth), '
         f'0 to {CLIP_SOLID} (SOLID-growth), growth rate ≥ {XL[0]:.2f} log10 doublings/yr; statistics use all data. '
         f'{MODE.get(CUR_MODE, CUR_MODE)}.')
fig.text(0.5, 0.005, wrap_caption(CAP13), ha='center', va='bottom', fontsize=6.3, color='dimgrey')
fig.subplots_adjust(left=0.115, right=0.955, top=0.90, bottom=0.20)
fig.savefig(OUT + "Fig13_twoaxis_density.png", bbox_inches='tight', dpi=300)
fig.savefig(OUT + "Fig13_twoaxis_density.pdf", bbox_inches='tight')
plt.close(fig)
log('saved Fig13_twoaxis_density.png / .pdf')


# ======================================================================
# Fig 14: aggressiveness gradient within baseline strata (approach A)
# ======================================================================
df['baseline'] = np.where(df['btype'] == 0, 'Baseline pGGN', 'Baseline PSN')
STRATA = ['Baseline pGGN', 'Baseline PSN']
STRAT_COL = {'Baseline pGGN': '#7FB3D5', 'Baseline PSN': '#E59866'}

sub = df.dropna(subset=['pn'])
log('')
log('=' * 70)
log('[Fig14] aggressiveness gradient within baseline strata (approach A)')
rows_jt = []
for st in STRATA:
    d = sub[sub['baseline'] == st]
    gs = [d.loc[d['group'] == g, 'pn'].values for g in GROUPS]
    z, p, ns = jt_test(gs)
    rows_jt.append((st, z, p, ns))
    log(f'  {st}: n={ns}  IAC%=' +
        ' / '.join(f'{100*d.loc[d["group"]==g,"IAC"].mean():.1f}' for g in GROUPS) +
        f'   JT z={z:+.2f}  {pfmt(p)}')
adj = bh([r[2] for r in rows_jt])
log('  BH-adjusted P: ' + ', '.join(f'{rows_jt[i][0]} {pfmt(adj[i])}' for i in range(len(rows_jt))))

fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), gridspec_kw={'width_ratios': [1.15, 1]})
ax = axes[0]
w = 0.36
for si, st in enumerate(STRATA):
    d = sub[sub['baseline'] == st]
    ks = [int(d.loc[d['group'] == g, 'IAC'].sum()) for g in GROUPS]
    ns = [int(d.loc[d['group'] == g, 'IAC'].notna().sum()) for g in GROUPS]
    ps = [k / n for k, n in zip(ks, ns)]
    lo, hi = zip(*[wilson(k, n) for k, n in zip(ks, ns)])
    x = np.arange(3) + (si - 0.5) * w
    ax.bar(x, ps, width=w, color=STRAT_COL[st], edgecolor='black', linewidth=0.6, label=st)
    ax.errorbar(x, ps, yerr=[np.array(ps) - np.array(lo), np.array(hi) - np.array(ps)],
                fmt='none', ecolor='black', elinewidth=0.7, capsize=2)
    for xi, p_, k, n, hi_ in zip(x, ps, ks, ns, hi):
        # Anchored to the **top of the error bar**, not the bar top: with the old
        # p_+0.045, a long error bar (the two SOLID-growth bars reach 0.95) let the
        # vertical line run through the middle of the text.
        # Same convention as Fig16 (ED_X3).
        ax.text(xi, max(p_, hi_) + 0.040, f'{100*p_:.1f}%', ha='center', fontsize=7,
                fontweight='bold')
        ax.text(xi, -0.055, f'{k}/{n}', ha='center', fontsize=6.3, color='dimgrey')
    yy = 1.13 if si == 1 else 1.05
    ax.text(1.0, yy, f'Jonckheere–Terpstra z = {rows_jt[si][1]:+.2f}, {pfmt(adj[si])}',
            ha='center', fontsize=7.4, color=STRAT_COL[st], fontweight='bold')
ax.set_xticks(np.arange(3))
ax.set_xticklabels(GROUPS, fontsize=8)
ax.set_ylim(-0.09, 1.22)
ax.set_yticks(np.arange(0, 1.01, 0.25))
ax.set_yticklabels([f'{int(v*100)}' for v in np.arange(0, 1.01, 0.25)])
ax.set_ylabel('Invasive adenocarcinoma (%)', fontsize=8.5)
ax.legend(frameon=False, fontsize=7.2, loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=2)
ax.text(-0.20, 1.08, 'A', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')
ax.set_title('Within every baseline density stratum', fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax = axes[1]
d = sub[sub['baseline'] == 'Baseline pGGN']
bottom = np.zeros(3)
for gi, (gl, gc) in enumerate(zip(GRADE_LAB, GRADE_COL)):
    vals = np.array([100 * (d.loc[d['group'] == g, 'pn'] == gi).mean() for g in GROUPS])
    # A stratum below 4% cannot fit its label inside the bar, so its actual n goes into the legend and the reader can still see that it exists
    _cnt = int((d['pn'] == gi).sum())
    ax.bar(np.arange(3), vals, bottom=bottom, color=gc, edgecolor='white',
           linewidth=0.6, width=0.62,
           label=(f'{gl} (n = {_cnt})' if 0 < max(vals) < 4 else gl))
    for xi, (v, b) in enumerate(zip(vals, bottom)):
        if v >= 4:
            ax.text(xi, b + v / 2, f'{v:.1f}', ha='center', va='center', fontsize=6.6,
                    color='white' if gi == 4 else 'black')
    bottom += vals
for xi, g in enumerate(GROUPS):
    ax.text(xi, 101.5, f'n = {int((d["group"] == g).sum())}', ha='center', fontsize=7, color='dimgrey')
ax.set_xticks(np.arange(3))
ax.set_xticklabels(GROUPS, fontsize=8)
ax.set_ylim(0, 108)
ax.set_ylabel('Resected nodules (%)', fontsize=8.5)
ax.legend(frameon=False, fontsize=7, ncol=5, loc='lower center', bbox_to_anchor=(0.5, -0.29),
          columnspacing=0.9, handlelength=1.1)
ax.text(-0.20, 1.08, 'B', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')
ax.set_title('Baseline pGGN only: full invasiveness spectrum', fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

CAP14 = ('Bars in A: invasive adenocarcinoma (IAC) rate with Wilson 95% confidence intervals; counts '
         '(IAC / all resected) below each bar; P values Benjamini–Hochberg-adjusted across the two strata. '
         'Baseline density was read on the first CT (pGGN, pure ground-glass nodule; PSN, part-solid nodule, '
         'with the heterogeneous ground-glass category merged into PSN). '
         'B: pathology grade distribution of the resected baseline-pGGN nodules; a grade occupying < 4% of a '
         'bar is not labelled inside the bar, and its nodule count is given in the legend instead.')
fig.text(0.5, 0.005, wrap_caption(CAP14), ha='center', va='bottom', fontsize=6.3, color='dimgrey')
fig.tight_layout(rect=[0, 0.085, 1, 1])
fig.savefig(OUT + "Fig14_baseline_stratified.png", bbox_inches='tight', dpi=300)
fig.savefig(OUT + "Fig14_baseline_stratified.pdf", bbox_inches='tight')
plt.close(fig)
log('saved Fig14_baseline_stratified.png / .pdf')


# ======================================================================
# Fig 15: three-group pairwise comparison, unadjusted vs adjusted for
# baseline 3D CTR + baseline volume
#   three rows = the three pairwise contrasts; two dots = unadjusted (open grey)
#   and adjusted (filled red). The adjustment terms are continuous covariates on a
#   different scale from the group contrasts, so their estimates are written to the
#   results table rather than drawn as a row.
# ======================================================================
pool = df.dropna(subset=['pn', 'CTR0_3d', 'V0_3d']).copy()
pool = pool[pool['V0_3d'] > 0]
D = pd.DataFrame({
    'y_bin': (pool['pn'] == 4).astype(int).values,
    'y_ord': pool['pn'].astype(int).values,
    'grp': pool['group'].values,
    'CTR0': pool['CTR0_3d'].astype(float).values,
    'log2V0': np.log2(pool['V0_3d'].astype(float).values),
}, index=pool.index)

# three rows: the two groups versus Stable + SOLID vs GGO with GGO-growth as reference
ROWS = [('GGO-growth vs Stable', 'Stable', 'GGO-growth'),
        ('SOLID-growth vs Stable', 'Stable', 'SOLID-growth'),
        ('SOLID-growth vs GGO-growth', 'GGO-growth', 'SOLID-growth')]
log('')
log('=' * 70)
log(f'[Fig15] three-group pairwise comparison, unadjusted vs adjusted for baseline 3D CTR+volume  n = {len(D)}  '
    f'(group distribution {pool["group"].value_counts().to_dict()})')


def design(ref):
    """Build the dummy-variable matrix for the given reference group."""
    cols = [g for g in GROUPS if g != ref]
    X = pd.DataFrame({f'grp_{g}': (D['grp'] == g).astype(float) for g in cols})
    return cols, X


def extract(r, key):
    names = list(r.params.index)
    ci = np.asarray(r.conf_int())
    i = names.index(key)
    p = r.pvalues.iloc[i] if hasattr(r.pvalues, 'iloc') else r.pvalues[i]
    return np.exp(r.params.iloc[i]), np.exp(ci[i][0]), np.exp(ci[i][1]), float(p)


res = {}
for ref in GROUPS:
    cols, X = design(ref)
    for lab, extra in [('unadj', []), ('adj', ['CTR0', 'log2V0'])]:
        Xf = pd.concat([X, D[extra]], axis=1) if extra else X
        res[('bin', lab, ref)] = sm.Logit(D['y_bin'], sm.add_constant(Xf)).fit(disp=0)
        res[('ord', lab, ref)] = OrderedModel(D['y_ord'], Xf, distr='logit').fit(method='bfgs', disp=0)

for lab, ref, key in ROWS:
    for kind in ('bin', 'ord'):
        e0 = extract(res[(kind, 'unadj', ref)], f'grp_{key}')
        e1 = extract(res[(kind, 'adj', ref)], f'grp_{key}')
        log(f'  [{kind}] {lab:29s} unadj OR={e0[0]:7.2f} ({e0[1]:.2f}-{e0[2]:.2f}) P={e0[3]:.2g}'
            f'   |  adj OR={e1[0]:7.2f} ({e1[1]:.2f}-{e1[2]:.2f}) P={e1[3]:.2g}')
# Estimates of the adjustment terms themselves (no figure row; used by the results table and caption)
for kind, nm in [('bin', 'IAC'), ('ord', 'ordered grade')]:
    e = extract(res[(kind, 'adj', 'Stable')], 'CTR0')
    log(f'  [{kind}] Baseline 3D CTR per +0.10 : OR={e[0]:.2f} ({e[1]:.2f}-{e[2]:.2f}) P={e[3]:.2g}')
    e = extract(res[(kind, 'adj', 'Stable')], 'log2V0')
    log(f'  [{kind}] Baseline volume per doubling: OR={e[0]:.2f} ({e[1]:.2f}-{e[2]:.2f}) P={e[3]:.2g}')

# Forest-plot layout: the left axis area holds only the dots and CIs; the numeric columns sit to the right of the plot area, so long intervals cannot squeeze the text into the panel gap
XLIM = (0.10, 200.0)
XTXT = XLIM[1] * 1.06
fig, axes = plt.subplots(1, 2, figsize=(10.2, 2.9), sharey=True)
for ai, (kind, ttl) in enumerate([('bin', 'Outcome: invasive adenocarcinoma'),
                                  ('ord', 'Outcome: ordered invasiveness grade')]):
    ax = axes[ai]
    ys = np.arange(len(ROWS))[::-1]
    for k, (lab, ref, key) in enumerate(ROWS):
        y = ys[k]
        for lname, off, filled in [('unadj', +0.16, False), ('adj', -0.16, True)]:
            orv, lo, hi, pv = extract(res[(kind, lname, ref)], f'grp_{key}')
            col = '#B03A2E' if filled else '#7F8C8D'
            ax.plot([lo, hi], [y + off] * 2, lw=1.2, color=col, solid_capstyle='butt')
            ax.plot([orv], [y + off], marker='o', ms=5.0,
                    mfc=col if filled else 'white', mec=col, mew=1.1, zorder=5)
            ax.text(XTXT, y + off, f'{orv:.2f} ({lo:.2f}–{hi:.2f})   P = {pv:.2g}',
                    fontsize=6.4, va='center', ha='left', color=col, clip_on=False)
    ax.axvline(1.0, ls='--', lw=0.8, c='grey')
    ax.set_xscale('log')
    ax.set_xlim(*XLIM)
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200])
    ax.set_xticklabels(['0.1', '0.2', '0.5', '1', '2', '5', '10', '20', '50', '100', '200'],
                       fontsize=7.2)
    ax.set_yticks(ys)
    if ai == 0:
        ax.set_yticklabels([r[0] for r in ROWS], fontsize=8.0)
    else:
        ax.tick_params(axis='y', labelleft=False)   # with sharey, set_yticklabels([]) clears A as well
    ax.set_xlabel('Odds ratio (log scale; the two groups are each other’s reciprocal)', fontsize=8)
    ax.set_title(ttl, fontsize=9)
    ax.text(-0.185 if ai == 0 else -0.04, 1.14, 'AB'[ai], transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='bottom')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(-0.55, len(ROWS) - 0.30)

axes[0].legend(handles=[Line2D([], [], marker='o', ls='none', mfc='white', mec='#7F8C8D',
                               mew=1.1, ms=5, label='Unadjusted'),
                        Line2D([], [], marker='o', ls='none', mfc='#B03A2E', mec='#B03A2E',
                               ms=5, label='Adjusted for baseline 3-D CTR and baseline volume')],
               frameon=False, fontsize=6.9, ncol=2, loc='upper center',
               bbox_to_anchor=(0.5, -0.26))
CAP15 = (f'Each row is a pairwise comparison between trajectory groups; the same fit with the roles reversed '
         f'gives the reciprocal (for example, Stable vs GGO-growth = 1/3.20). Open circles, unadjusted; filled '
         f'circles, adjusted for baseline 3-D CT ratio and baseline volume. Panel A, logistic regression for '
         f'invasive adenocarcinoma; panel B, proportional-odds ordered logit for the ordered invasiveness grade '
         f'(benign < AAH < AIS < MIA < IAC). n = {len(D)} resected nodules with baseline 3-D CT metrics. '
         f'Adjustment left every comparison essentially unchanged. The two adjustment terms themselves are '
         f'reported separately rather than plotted here, because they are continuous covariates '
         f'whose units are not comparable with a between-group odds ratio.')
fig.text(0.5, -0.02, wrap_caption(CAP15), ha='center', va='top', fontsize=6.3, color='dimgrey')
fig.subplots_adjust(left=0.155, right=0.775, top=0.84, bottom=0.26, wspace=0.60)
fig.savefig(OUT + "Fig15_ctr_adjusted_forest.png", bbox_inches='tight', dpi=300)
fig.savefig(OUT + "Fig15_ctr_adjusted_forest.pdf", bbox_inches='tight')
plt.close(fig)
log('saved Fig15_ctr_adjusted_forest.png / .pdf')


# ======================================================================
# Fig 16: aggressiveness gradient within last-CT strata
#   Definition caveat: the last-CT state is a downstream consequence of the
#   trajectory (a mediator), not a pre-existing factor (a confounder). So this is
#   not confounding adjustment but a specificity test: does the trajectory carry
#   information beyond the "last imaging reading"? In the last-CT pGGN stratum
#   SOLID-growth has only 1 nodule (SOLID-growth densifies by definition), so it
#   cannot be estimated and is labelled as such.
# ======================================================================
FIN_COL = {'Final pGGN': '#9DC3E6', 'Final PSN': '#E59866'}
fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), gridspec_kw={'width_ratios': [1.15, 1]})

# --- A: last-CT PSN stratum, three groups ---
ax = axes[0]
sA = sub[sub['ftype'] != 0]
ks = [int(((sA['group'] == g) & (sA['pn'] == 4)).sum()) for g in GROUPS]
ns = [int((sA['group'] == g).sum()) for g in GROUPS]
ps = [k / n if n else np.nan for k, n in zip(ks, ns)]
lo, hi = zip(*[wilson(k, n) for k, n in zip(ks, ns)])
pos = np.arange(3)
ax.bar(pos, ps, width=0.56, color=[COLORS[g] for g in GROUPS], edgecolor='black', linewidth=0.6)
ax.errorbar(pos, ps, yerr=[np.array(ps) - np.array(lo), np.array(hi) - np.array(ps)],
            fmt='none', ecolor='black', elinewidth=0.7, capsize=2)
for xi, p_, k, n, hi_ in zip(pos, ps, ks, ns, hi):
    ax.text(xi, max(p_, hi_) + 0.040, f'{100*p_:.1f}%', ha='center', fontsize=7.2, fontweight='bold')
    ax.text(xi, -0.055, f'{k}/{n}', ha='center', fontsize=6.3, color='dimgrey')
zA, pA, _ = jt_test([sA.loc[sA['group'] == g, 'pn'].values for g in GROUPS])
ax.text(0.5, 0.955, f'Jonckheere–Terpstra z = {zA:+.2f}, {pfmt(pA)}', transform=ax.transAxes,
        ha='center', va='top', fontsize=7.6, color='#8A4B1E', fontweight='bold')
ax.set_xticks(pos)
ax.set_xticklabels(GROUPS, fontsize=8)
ax.set_ylim(-0.09, 1.20)
ax.set_yticks(np.arange(0, 1.01, 0.25))
ax.set_yticklabels([f'{int(v*100)}' for v in np.arange(0, 1.01, 0.25)])
ax.set_ylabel('Invasive adenocarcinoma (%)', fontsize=8.5)
ax.text(-0.20, 1.08, 'A', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')
ax.set_title('Nodules still part-solid at the last CT', fontsize=9, pad=6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
log('')
log(f'[Fig16] last-CT PSN stratum n={ns}  IAC%={[round(100*p_,1) for p_ in ps]}  JT z={zA:+.2f} P={pA:.2g}')

# --- B: last-CT pGGN stratum, only Stable and GGO-growth are estimable ---
ax = axes[1]
sB = sub[sub['ftype'] == 0]
GB = ['Stable', 'GGO-growth']
ksB = [int(((sB['group'] == g) & (sB['pn'] == 4)).sum()) for g in GB]
nsB = [int((sB['group'] == g).sum()) for g in GB]
psB = [k / n if n else np.nan for k, n in zip(ksB, nsB)]
loB, hiB = zip(*[wilson(k, n) for k, n in zip(ksB, nsB)])
posB = np.arange(3)
ax.bar(posB[:2], psB, width=0.56, color=[COLORS[g] for g in GB], edgecolor='black', linewidth=0.6)
ax.errorbar(posB[:2], psB, yerr=[np.array(psB) - np.array(loB), np.array(hiB) - np.array(psB)],
            fmt='none', ecolor='black', elinewidth=0.7, capsize=2)
for xi, p_, k, n, hi_ in zip(posB[:2], psB, ksB, nsB, hiB):
    ax.text(xi, max(p_, hi_) + 0.040, f'{100*p_:.1f}%', ha='center', fontsize=7.2, fontweight='bold')
    ax.text(xi, -0.055, f'{k}/{n}', ha='center', fontsize=6.3, color='dimgrey')
nS = int((sB['group'] == 'SOLID-growth').sum())
ax.text(posB[2], 0.5, f'SOLID-growth\nn = {nS}\n(not estimable)', ha='center', va='center',
        fontsize=6.8, color='dimgrey')
zB, pB, _ = jt_test([sB.loc[sB['group'] == g, 'pn'].values for g in GB])
ax.text(0.42, 0.955, f'Jonckheere–Terpstra z = {zB:+.2f}, {pfmt(pB)}', transform=ax.transAxes,
        ha='center', va='top', fontsize=7.6, color='#2E75B6', fontweight='bold')
ax.set_xticks(posB)
ax.set_xticklabels(GROUPS, fontsize=8)
ax.set_ylim(-0.09, 1.20)
ax.set_yticks(np.arange(0, 1.01, 0.25))
ax.set_yticklabels([f'{int(v*100)}' for v in np.arange(0, 1.01, 0.25)])
ax.set_ylabel('Invasive adenocarcinoma (%)', fontsize=8.5)
ax.text(-0.20, 1.08, 'B', transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')
ax.set_title('Nodules still pure ground-glass at the last CT', fontsize=9, pad=6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
_pct_solid_psn = 100 * (sub.loc[sub['group'] == 'SOLID-growth', 'ftype'] != 0).mean()
log(f'[Fig16] last-CT pGGN stratum n={nsB}(Stable/GGO) + SOLID n={nS}  IAC%={[round(100*p_,1) for p_ in psB]}'
    f'  JT z={zB:+.2f} P={pB:.2g}')
log(f'[Fig16] share of SOLID-growth nodules still pGGN at the last CT = {100-_pct_solid_psn:.1f}%')

CAP16 = ('Density at the last CT was used to split the same resected nodules into two strata. This is a '
         'specificity check, not a confounding control: density at the last CT lies downstream of the '
         f'trajectory, so conditioning on it can only attenuate the gradient, never manufacture one. '
         f'A: nodules that were part-solid at the last CT (n = {ns[0]}/{ns[1]}/{ns[2]} for Stable / GGO-growth / '
         f'SOLID-growth). B: nodules that were still pure ground-glass at the last CT (n = {nsB[0]}/{nsB[1]}); '
         f'only {nS} SOLID-growth nodule fell in this stratum — as expected, because densification is how that '
         f'trajectory is defined ({_pct_solid_psn:.1f}% of SOLID-growth nodules were part-solid at the last CT) — '
         f'so the comparison is restricted to Stable vs GGO-growth. Bars: invasive adenocarcinoma rate with '
         f'Wilson 95% confidence intervals; counts (IAC / all resected) below each bar.')
fig.text(0.5, -0.02, wrap_caption(CAP16), ha='center', va='top', fontsize=6.3, color='dimgrey')
fig.tight_layout()
fig.savefig(OUT + "Fig16_finalCT_stratified.png", bbox_inches='tight', dpi=300)
fig.savefig(OUT + "Fig16_finalCT_stratified.pdf", bbox_inches='tight')
plt.close(fig)
log('saved Fig16_finalCT_stratified.png / .pdf')

