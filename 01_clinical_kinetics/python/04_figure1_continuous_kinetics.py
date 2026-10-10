#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Continuous kinetics across the three growth patterns (five panels)
#
# Purpose   : Redraw the continuous-kinetics panel at the final size: VDT, MDT of
#             total mass, MDT of solid mass, annualized CTR and follow-up years,
#             in the style of the corresponding panel of 03_figure1_panels.py.
#
# Inputs    : ./kinetics.csv (built by 01_build_nodule_kinetics_table.py)
#
# Outputs   : ./figure1_panels/Fig1f_continuous_kinetics.png / .pdf
#             ./figure1_panels/Fig1f_statistics.csv (read back by 03)
#             ./04_figure1_continuous_kinetics.log
#
# Usage     : PROJECT_ROOT=/path/to/data python3 04_figure1_continuous_kinetics.py
#             Normally launched by 03_figure1_panels.py, which reads back
#             figure1_panels/Fig1f_statistics.csv; it can also run on its own.
# -----------------------------------------------------------------------------

import os, csv
import numpy as np
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import textwrap
from matplotlib import rcParams
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__)) + "/"
OUT = os.path.dirname(os.path.abspath(__file__)) + "/figure1_panels/"
os.makedirs(OUT, exist_ok=True)

GROUPS = ['Stable', 'GGO-growth', 'SOLID-growth']
COLORS = {'Stable': '#4E79A7', 'GGO-growth': '#F28E2B', 'SOLID-growth': '#E15759'}

rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'],
                 'font.size': 8.5, 'axes.linewidth': 0.7, 'xtick.direction': 'out',
                 'ytick.direction': 'out', 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                 'figure.dpi': 300, 'savefig.bbox': 'tight'})

LOG, STATS = [], []


def log(s=''):
    print(s); LOG.append(str(s))


def num(x):
    try:
        s = str(x).strip()
        return None if s == '' else float(s)
    except (TypeError, ValueError):
        return None


def clean_group(g):
    s = str(g).strip()
    return 'GGO-growth' if s == 'GK253GO-growth' else s


def bh(ps):
    ps = np.asarray(ps, float); k = len(ps)
    order = np.argsort(ps); adj = np.empty(k)
    adj[order] = np.minimum.accumulate((ps[order] * k / np.arange(1, k + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)


def star(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


# ======================================================================
# Read the data
# ======================================================================
K = list(csv.DictReader(open(BASE + "kinetics.csv", encoding="utf-8-sig")))
for r in K:
    r['group'] = clean_group(r['group'])
K = [r for r in K if r['group'] in GROUPS]
N_TOT = Counter(r['group'] for r in K)

log("=" * 78)
log("[f] Continuous kinetics across the three groups - five panels (VDT / MDT total / MDT solid / annualized CTR / follow-up)")
log("=" * 78)
log("  Nodules per group: " + "  ".join(f"{g}={N_TOT[g]}" for g in GROUPS))

# (column name, y-axis label, log10?, 99th-percentile clip?, display window)
METRICS = [
    ('VDT',         'VDT (days, log$_{10}$)',              True,  True,  None),
    ('MDT3d_total', 'MDT, total mass (days, log$_{10}$)',  True,  True,  None),
    ('MDT3d_solid', 'MDT, solid mass (days, log$_{10}$)',  True,  True,  None),
    ('ann_CTR3d',   'Annualized CTR change (/yr)',         False, False,
     {'Stable': (-0.20, 0.20), 'GGO-growth': (-0.20, 0.20), 'SOLID-growth': (0.00, 0.30)}),
    ('years',       'Follow-up (years)',                   False, False, None),
]
CAP_W = 160     # caption wrap width (characters). See the fig.text call below.
LOGM = ('VDT', 'MDT3d_total', 'MDT3d_solid')          # only > 0 values take a logarithm

# ---- Follow-up panel inclusion floor (fifth cell only) ----------------------
# The Stable group is restricted to follow-up >= 1 year: a nodule followed for
# 0.04 years has no right to be called "Stable" - with such a short observation
# window, "did not grow" does not imply "stable". The two growth groups are not
# restricted: they already grew, so a short follow-up still holds, and dropping
# them would only throw away true signal.
# The first four cells keep the full data: their definition is a rate / first-to-
# last change, independent of the length of the observation window.
# The restriction applies to both full and disp - Kruskal-Wallis and the pairwise
# tests are computed on the plotted set, consistent with the existing convention
# that "shown n = analyzable n".
MIN_VAL = {'years': {'Stable': 1.0}}

# ---- Values: full = complete data (for statistics + KW), data = actually plotted --
FULL, DISP = [], []
for m, _, _, _, win in METRICS:
    full, disp = [], []
    gmin = MIN_VAL.get(m, {})
    for g in GROUPS:
        vals = [num(r.get(m)) for r in K if r['group'] == g]
        vals = np.array([v for v in vals if v is not None], float)
        if m in LOGM:
            vals = vals[vals > 0]
        if g in gmin:
            vals = vals[vals >= gmin[g]]
        full.append(vals)
        d = np.log10(vals) if m in LOGM else vals.copy()
        if m in LOGM:
            d = d[d <= np.percentile(d, 99)]
        if win is not None:
            lo_w, hi_w = win[g]
            d = d[(d >= lo_w) & (d <= hi_w)]
        disp.append(np.array(d, float))
    FULL.append(full); DISP.append(disp)

# ---- Pairwise Mann-Whitney + BH (**complete data**) --------------------------
pair_ps = []
for full in FULL:
    ps = []
    for a, b in ((0, 1), (0, 2), (1, 2)):
        ps.append(stats.mannwhitneyu(full[a], full[b], alternative='two-sided')[1])
    pair_ps.append(bh(ps))

# ======================================================================
# Canvas: three cells in the top row, two half-width cells in the bottom row
# ======================================================================
# Margins tightened from the gridspec defaults 0.125/0.9 to 0.055/0.998: without
# this, each top-row cell is only 1.4 in wide and "GGO-growth" and "SOLID-growth"
# merge into one word.
#
# The five cells are **equal width**: a 6-column grid, each cell spanning 2
# columns. The previous version gave the two bottom cells 3 columns each
# (3.36 in), 60% wider than the top row (2.09 in), which stretched the violins
# horizontally and made the five-panel figure look like two different layouts.
# Now the bottom cells span 2 columns each, offset by one column (1:3 and 3:5)
# and centred, so all five cells are 1.98 in and close to square.
# Canvas width 7.6 -> 7.2 -> 6.6 in. The lower bound is set by the x-axis labels:
# after the three columns are evenly split, each cell slot must fit the widest
# label "SOLID-growth" (measured 39.2 pt at 6.2 pt). At 6.6 in the slot is
# 43.2 pt, 4.0 pt of margin; at 6.4 in only 2.7 pt remain and the two labels would
# touch. A narrower canvas would require rotating the group names by 30 deg.
fig = plt.figure(figsize=(6.6, 5.4))
gs = fig.add_gridspec(2, 6, left=0.062, right=0.998, top=0.93, bottom=0.10,
                      hspace=0.52, wspace=0.55)
axes = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]),
        fig.add_subplot(gs[1, 1:3]), fig.add_subplot(gs[1, 3:5])]

for pi, (m, ylab, logscale, clip, win) in enumerate(METRICS):
    ax, data, full = axes[pi], DISP[pi], FULL[pi]
    # Must be re-fetched inside this loop: after the previous loop gmin is left
    # on the last item (years); reusing it would wrongly report the "Stable>=1"
    # floor in the log of the first four cells.
    gmin = MIN_VAL.get(m, {})
    pos = [0, 1, 2]
    for k, g in enumerate(GROUPS):
        vp = ax.violinplot([data[k]], positions=[pos[k]], widths=0.62,
                           showextrema=False, showmedians=False)
        for pc in vp['bodies']:
            pc.set_facecolor(COLORS[g]); pc.set_alpha(0.45)
            pc.set_edgecolor('black'); pc.set_linewidth(0.6)
    bp = ax.boxplot(data, positions=pos, widths=0.13, patch_artist=True,
                    showfliers=False, medianprops=dict(color='black', lw=1.2))
    for patch in bp['boxes']:
        patch.set_facecolor('white'); patch.set_alpha(0.9); patch.set_linewidth(0.7)

    # n uses "analyzable / group total".
    # The 99th-percentile clipping of the log panels and the CTR display window
    # only affect the plotted range, not the analyzable n, so the labels do not
    # change with them; the CTR window is noted separately in the figure footer.
    # Exception: the follow-up panel's Stable group has a >= 1 year floor that
    # **does** change the analyzable n, so that cell shows 514/575 and the reason
    # is noted in the figure footer.
    ax.set_xticks(pos)
    ax.set_xticklabels([f"{g}\n(n={len(full[k])}/{N_TOT[g]})" for k, g in enumerate(GROUPS)],
                       fontsize=6.2, linespacing=1.4)
    ax.set_ylabel(ylab, fontsize=7.6)

    _, pkw = stats.kruskal(*full)                      # KW on the complete data
    ax.text(0.5, 1.045, f'KW P = {pkw:.1e}', transform=ax.transAxes, ha='center',
            va='bottom', fontsize=6.2, color='dimgrey')

    # The order of the names here must match the order in which pair_ps was built
    # above, ((0,1),(0,2),(1,2)); the assertion below guards the correspondence.
    p_01, p_02, p_12 = pair_ps[pi]
    assert np.allclose([p_01, p_02, p_12], pair_ps[pi]), "bracket p values do not match the order of the comparison pairs"
    top = max(np.max(d) for d in data if len(d))
    lo = min(np.min(d) for d in data if len(d))
    h = (top - lo) * 0.035
    for (x1, x2), p in (((0, 1), p_01), ((0, 2), p_02), ((1, 2), p_12)):
        y = top + (top - lo) * (0.08 if (x1, x2) != (0, 2) else 0.21)
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.7, color='#333333')
        ax.text((x1 + x2) / 2, y + h, star(p), ha='center', va='bottom', fontsize=7.0)
    ax.set_ylim(lo - (top - lo) * 0.08, top + (top - lo) * 0.34)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

    cut = "" if win is None else f"  {sum(len(f) - len(d) for f, d in zip(full, data))} nodules outside the display window (statistics unaffected)"
    if gmin:
        cut += "  inclusion floor " + " ".join(
            f"{g}>={v:g}" for g, v in gmin.items()) + " (already counted in analyzable n)"
    log(f"  {m:14s} KW P={pkw:.2e}  shown_n=" + "/".join(str(len(d)) for d in data)
        + "  analyzable_n=" + "/".join(str(len(f)) for f in full) + cut)

    STATS.append(dict(panel='Fig1f', metric=m, ylab=ylab,
                      kw_p=f"{pkw:.6e}",
                      bh_p_01=f"{p_01:.6g}", bh_p_02=f"{p_02:.6g}", bh_p_12=f"{p_12:.6g}",
                      n_shown="/".join(str(len(d)) for d in data),
                      n_analyzable="/".join(str(len(f)) for f in full),
                      n_group_total="/".join(str(N_TOT[g]) for g in GROUPS),
                      n_outside_window=sum(len(f) - len(d) for f, d in zip(full, data))))

# The panel letter must use **figure coordinates**, not axes fractions. Pinned at
# x=-0.34 on axes[0], the left axis edge is at 0.062*6.6=0.41 in and
# -0.34*1.80 = -0.61 in pushes the letter to -0.20 in, outside the canvas;
# savefig.bbox='tight' then extends the whole canvas 0.2 in to the left - the
# panels are squeezed but the canvas does not follow.
fig.text(0.010, 0.978, 'f', fontsize=12, fontweight='bold', va='bottom')

# The display window must be stated in the figure: the annualized-CTR panel plots
# only the 535/268/49 nodules inside the window; without saying so, a reader would
# assume that is the whole distribution of the three groups
n_cut = sum(len(f) - len(d) for f, d in zip(FULL[3], DISP[3]))
# The caption must wrap inside the canvas. savefig.bbox='tight' sizes the canvas
# by the **widest content**: these two sentences total about 400 characters,
# about 7.8 in wide at 5.6 pt, wider than the 6.6 in canvas - without wrapping,
# shrinking the panels to 1.80 in does not help, the figure is still 7.8 in wide.
_cap = (f"Annualized CTR panel: {n_cut} nodules outside the display window "
        f"(|ΔCTR/yr| > 0.20 for Stable and GGO-growth; < 0 or > 0.30 for SOLID-growth) "
        f"are not shown. All statistics were computed on complete data. "
        f"Follow-up panel: the Stable group is restricted to nodules followed for "
        f"at least 1 year ({N_TOT['Stable'] - len(FULL[4][0])} of {N_TOT['Stable']} "
        f"excluded); a nodule observed for less than a year cannot be called stable. "
        f"The two growth groups are unconstrained.")
fig.text(0.010, -0.012, textwrap.fill(_cap, width=CAP_W),
         fontsize=5.6, color='dimgrey', ha='left', va='top', linespacing=1.5)

fig.savefig(OUT + "Fig1f_continuous_kinetics.png", dpi=300)
fig.savefig(OUT + "Fig1f_continuous_kinetics.pdf")
plt.close(fig)
log(f"  saved Fig1f_continuous_kinetics.png / .pdf")

with open(OUT + "Fig1f_statistics.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(STATS[0])); w.writeheader(); w.writerows(STATS)
log("  statistics -> figure1_panels/Fig1f_statistics.csv")

with open(BASE + "04_figure1_continuous_kinetics.log", "w", encoding="utf-8") as f:
    f.write("\n".join(LOG) + "\n")
