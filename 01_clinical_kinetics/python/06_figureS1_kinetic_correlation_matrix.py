#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Correlation structure among the kinetic metrics
#
# Purpose   : Pairwise Spearman correlations among the four kinetic metrics (95% CI
#             and R^2), the variance decomposition of mass change into volume and
#             density, and a robustness check across the three VDT definitions;
#             draws the correlation heat map and scatter panels.
#
# Inputs    : ./kinetics.csv (built by 01_build_nodule_kinetics_table.py); for the
#             robustness table, sibling folders *VDT_MDT_continuous_kinetics_*
#
# Outputs   : ./SuppTable_kinetic_correlation.csv, ./Fig17_kinetic_correlation.png/.pdf,
#             ./kinetic_correlation_report.md
#
# Usage     : PROJECT_ROOT=/path/to/data python3 06_figureS1_kinetic_correlation_matrix.py
# -----------------------------------------------------------------------------
"""
Correlation structure among the kinetic quantities + variance decomposition of
mass = volume x density
Output: kinetic_correlation_report.md / SuppTable_kinetic_correlation.csv /
        Fig17_kinetic_correlation.png + .pdf
Definition: this folder's kinetics.csv (MODE follows the rows)
"""
# -----------------------------------------------------------------------------

import os, math
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = HERE + os.sep

D = pd.read_csv(os.path.join(HERE, 'kinetics.csv'))
for c in ['VDT','VDT2d','VDT3d','MDT3d_total','MDT3d_solid','ann_CTR3d',
          'V0_3d','V1_3d','MT0','MT1','MS0','MS1','CTR0_3d','CTR1_3d']:
    D[c] = pd.to_numeric(D[c], errors='coerce')

GROUPS = ['Stable', 'GGO-growth', 'SOLID-growth']
COLORS = {'Stable': '#4C72B0', 'GGO-growth': '#55A868', 'SOLID-growth': '#C44E52'}
METRICS = [('VDT','VDT','Volume doubling time'),
           ('MDT3d_total','MDT_total','Mass doubling time (total)'),
           ('MDT3d_solid','MDT_solid','Mass doubling time (solid)'),
           ('ann_CTR3d','ann_CTR','Annualized CTR change')]

def fisher_ci(r, n, alpha=0.05):
    if n < 5 or not np.isfinite(r) or abs(r) >= 1: return (np.nan, np.nan)
    z = np.arctanh(r); se = 1.0/np.sqrt(n-3); zc = stats.norm.ppf(1-alpha/2)
    return float(np.tanh(z-zc*se)), float(np.tanh(z+zc*se))

L = []
def p(s=''):
    L.append(s); print(s)

p("="*100)
p("Correlation structure of the kinetic metrics")
p("="*100)
p(f"data: {os.path.basename(HERE)}/kinetics.csv   total rows {len(D)}")
p("")

# ---------- Table 1: pairwise Spearman ----------
p("="*100); p("[Table 1] pairwise Spearman rank correlations among the four kinetic metrics (pairwise complete, with 95% CI and R2)"); p("="*100)
hdr = f"  {'':34s}" + ''.join(f"{lab:>22s}" for _, lab, _ in METRICS)
p(hdr)
tab_rows = []
for a, la, lf in METRICS:
    line = f"  {lf:34s}"
    for b, lb, lf2 in METRICS:
        if a == b:
            line += f"{'—':>22s}"; continue
        s = D[[a, b]].dropna()
        r = stats.spearmanr(s[a], s[b]).statistic
        line += f"{r:>+13.3f} (n={len(s):3d})"
        if a < b or True:
            lo, hi = fisher_ci(r, len(s))
            tab_rows.append(dict(metric_1=lf, metric_2=lf2, n=len(s),
                                 spearman_rho=round(r,3), ci_low=round(lo,3), ci_high=round(hi,3),
                                 r_squared=round(r*r,3),
                                 p_value=stats.spearmanr(s[a], s[b]).pvalue))
    p(line)
p("")
p("  Main reading (three key pairs):")
KEY = [('VDT','ann_CTR3d','VDT  x  annualized CTR'),
       ('MDT3d_total','ann_CTR3d','MDT total x  annualized CTR'),
       ('MDT3d_solid','ann_CTR3d','MDT solid x  annualized CTR')]
for a, b, lab in KEY:
    s = D[[a, b]].dropna(); r = stats.spearmanr(s[a], s[b]).statistic
    lo, hi = fisher_ci(r, len(s))
    p(f"    {lab:22s} rho = {r:+.3f}  95% CI [{lo:+.3f}, {hi:+.3f}]  n = {len(s):3d}  R2 = {r*r:.3f}")
p("    -> VDT and the densification axis are almost orthogonal; the more solid the MDT, the more it is 'pulled away' by the densification axis.")
p("")

# ---------- Table 2: variance decomposition of mass = volume x density ----------
p("="*100); p("[Table 2] Variance decomposition of mass change: mass = volume x density (an identity)"); p("="*100)
s = D[['V0_3d','V1_3d','MT0','MT1','ann_CTR3d']].dropna()
s = s[(s['V0_3d'] > 0) & (s['V1_3d'] > 0) & (s['MT0'] > 0) & (s['MT1'] > 0)]
dlV = np.log(s['V1_3d']/s['V0_3d'])
dlM = np.log(s['MT1']/s['MT0'])
dlR = np.log((s['MT1']/s['V1_3d'])/(s['MT0']/s['V0_3d']))    # density change
p(f"  fully decomposable n = {len(s)}")
p(f"  identity max|dlnM - (dlnV + dlnRho)| = {np.max(np.abs(dlM-(dlV+dlR))):.3e}   (numerical check)")
vV, vR, cVR = dlV.var(ddof=1), dlR.var(ddof=1), np.cov(dlV, dlR)[0,1]
vM = dlM.var(ddof=1)
p(f"  Var(dln Mass) = {vM:.5f}")
p(f"      = Var(dln Volume) {vV:.5f}  +  Var(dln Density) {vR:.5f}  +  2*Cov {2*cVR:.5f}")
p(f"  shares: volume {vV/vM*100:5.1f}%    density {vR/vM*100:5.1f}%    covariance {2*cVR/vM*100:5.1f}%")
p(f"  Spearman rho(dlnVolume, dlnDensity) = {stats.spearmanr(dlV, dlR).statistic:+.3f}")
p(f"  -> about {vR/vM*100:.0f}% of the variation in mass change comes from density change, which VDT cannot see at all.")
p("")

# ---------- Table 3: robustness across three definitions ----------
p("="*100); p("[Table 3] VDT x MDT agreement across the three VDT definitions (robustness check)"); p("="*100)
p(f"  {'definition':14s}{'VDT x MDT total':>20s}{'VDT x MDT solid':>20s}{'n':>10s}")
for f in ['VDT_MDT_continuous_kinetics_VDT2d','VDT_MDT_continuous_kinetics_VDT2d_fill3d','VDT_MDT_continuous_kinetics_VDT3d']:
    fp = os.path.join(os.path.dirname(HERE), f, 'kinetics.csv')
    if not os.path.exists(fp): continue
    d2 = pd.read_csv(fp)
    for c in ['VDT','MDT3d_total','MDT3d_solid']: d2[c] = pd.to_numeric(d2[c], errors='coerce')
    a = d2[['VDT','MDT3d_total']].dropna(); b = d2[['VDT','MDT3d_solid']].dropna()
    p(f"  {f.replace('VDT_MDT_continuous_kinetics_',''):14s}"
      f"{stats.spearmanr(a['VDT'],a['MDT3d_total']).statistic:>+13.3f}(n={len(a):3d})"
      f"{stats.spearmanr(b['VDT'],b['MDT3d_solid']).statistic:>+13.3f}(n={len(b):3d})")
p("")

# ---------- Export tables ----------
pd.DataFrame(tab_rows).to_csv(os.path.join(OUT, 'SuppTable_kinetic_correlation.csv'),
                              index=False, encoding='utf-8-sig')
p(f"  exported: SuppTable_kinetic_correlation.csv ({len(tab_rows)} pairs)")

# ================= Figures =================
plt.rcParams.update({'font.family':'Arial','font.size':7.2,'axes.linewidth':0.7,
                     'xtick.major.width':0.7,'ytick.major.width':0.7})
fig = plt.figure(figsize=(7.4, 4.9))
gs  = fig.add_gridspec(2, 3, hspace=0.52, wspace=0.40, left=0.085, right=0.985,
                       top=0.885, bottom=0.235)

# A: correlation heat map
axA = fig.add_subplot(gs[0, 0])
short = ['VDT', 'MDT\ntotal', 'MDT\nsolid', 'ann\nCTR']
M = np.full((4, 4), np.nan)
for i, (a, _, _) in enumerate(METRICS):
    for j, (b, _, _) in enumerate(METRICS):
        if i == j: continue
        ss = D[[a, b]].dropna()
        if len(ss) > 4: M[i, j] = stats.spearmanr(ss[a], ss[b]).statistic
im = axA.imshow(M, cmap='RdBu_r', vmin=-1, vmax=1)
axA.set_xticks(range(4)); axA.set_xticklabels(short, fontsize=6.4)
axA.set_yticks(range(4)); axA.set_yticklabels(short, fontsize=6.4)
for i in range(4):
    for j in range(4):
        if i == j:
            axA.text(j, i, '—', ha='center', va='center', fontsize=7, color='grey'); continue
        axA.text(j, i, f'{M[i,j]:+.2f}', ha='center', va='center', fontsize=6.6,
                 color='white' if abs(M[i, j]) > 0.55 else 'black')
axA.set_title('Pairwise Spearman $\\rho$', fontsize=7.6, pad=4)
axA.tick_params(length=2)
cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.04); cb.ax.tick_params(labelsize=6)
axA.text(-0.42, 1.22, 'a', transform=axA.transAxes, fontsize=10, fontweight='bold')

# B-D: three scatter panels
PANELS = [('VDT', 'ann_CTR3d', 'log10 VDT (days)', 'Annualized CTR change (/yr)', True),
          ('MDT3d_total', 'ann_CTR3d', 'log10 MDT total (days)', 'Annualized CTR change (/yr)', False),
          ('MDT3d_solid', 'ann_CTR3d', 'log10 MDT solid (days)', 'Annualized CTR change (/yr)', False)]
AXB = [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[1, 0])]
for k, ((a, b, xl, yl, logx), ax) in enumerate(zip(PANELS, AXB)):
    for g in GROUPS:
        sg = D[D['group'] == g][[a, b]].dropna()
        if not len(sg): continue
        x = np.log10(sg[a].clip(lower=1)) if logx else np.log10(sg[a].clip(lower=1))
        ax.scatter(x, sg[b], s=5, alpha=0.40, c=COLORS[g], edgecolors='none', label=g)
    ss = D[[a, b]].dropna()
    r = stats.spearmanr(ss[a], ss[b]).statistic
    ax.axhline(0, color='grey', lw=0.6, ls='--')
    ax.set_xlabel(xl, fontsize=6.9); ax.set_ylabel(yl, fontsize=6.9)
    ax.text(0.03, 0.965, f'$\\rho$ = {r:+.2f}   n = {len(ss)}', transform=ax.transAxes,
            ha='left', va='top', fontsize=6.8)
    ax.text(-0.30, 1.16, chr(98 + k), transform=ax.transAxes, fontsize=10, fontweight='bold')
    ax.tick_params(labelsize=6.4, length=2)
    if k == 0: ax.legend(fontsize=6, loc='upper right', frameon=False, markerscale=2.2,
                         handletextpad=0.3, borderpad=0.2)

# E: density-loading bars
axE = fig.add_subplot(gs[1, 1:])
labs = ['VDT\n(volume)', 'MDT total\n(volume$\\times$density)', 'MDT solid\n(solid vol$\\times$density)']
vals, ns = [], []
for a in ['VDT', 'MDT3d_total', 'MDT3d_solid']:
    ss = D[[a, 'ann_CTR3d']].dropna()
    vals.append(stats.spearmanr(ss[a], ss['ann_CTR3d']).statistic); ns.append(len(ss))
cols = ['#4C72B0', '#8E7CC3', '#C44E52']
bars = axE.bar(range(3), vals, color=cols, width=0.55, edgecolor='black', lw=0.6)
for i, (v, n) in enumerate(zip(vals, ns)):
    axE.text(i, v - 0.055 if v < 0 else v + 0.02, f'{v:+.2f}\n(n={n})', ha='center',
             va='top' if v < 0 else 'bottom', fontsize=6.6)
axE.axhline(0, color='black', lw=0.7)
axE.set_xticks(range(3)); axE.set_xticklabels(labs, fontsize=6.6)
axE.set_ylabel('Spearman $\\rho$ vs annualized CTR change', fontsize=6.9)
axE.set_ylim(-0.85, 0.30)
axE.tick_params(labelsize=6.4, length=2)
axE.text(-0.13, 1.14, 'e', transform=axE.transAxes, fontsize=10, fontweight='bold')
axE.set_title('Density loading of each growth metric', fontsize=7.6, pad=5)

# Was 0.175 - the x-axis title of panel d (gs[1,0], bottom=0.235),
# "log10 MDT solid (days)", had its baseline at ~0.187, below the top of the
# footnote, so the two lines overlapped. Moved down to 0.150: the footnote
# baseline is 0.031 and 0.15 in of the figure is still left at the bottom.
fig.text(0.5, 0.150,
    'a, Spearman rank correlations among the four time-normalized kinetic metrics (pairwise complete observations).\n'
    'b, Total-lesion volume doubling time against annualized CTR change. c, Total-mass doubling time against annualized\n'
    'CTR change. d, Solid-mass doubling time against annualized CTR change. e, Rank correlation of each growth metric with\n'
    'annualized CTR change, showing the progressive density loading of mass-based metrics. Mass is the product of volume\n'
    'and mean attenuation, so mass-based kinetics incorporate densification that volume-based kinetics do not.',
    ha='center', va='top', fontsize=5.9, color='dimgrey', linespacing=1.45)

fig.savefig(os.path.join(OUT, 'Fig17_kinetic_correlation.png'), dpi=600)
fig.savefig(os.path.join(OUT, 'Fig17_kinetic_correlation.pdf'))
p(f"  figure written: Fig17_kinetic_correlation.png / .pdf")

with open(os.path.join(OUT, 'kinetic_correlation_report.md'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(L))
print("  report: kinetic_correlation_report.md")
