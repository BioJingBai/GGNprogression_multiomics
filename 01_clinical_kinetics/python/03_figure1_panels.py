# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Clinical-kinetics panel assembly
#
# Purpose   : Build the clinical-kinetics panels from the analysis-ready kinetics
#             table: the Sankey transition, the baseline nodule-type composition,
#             the continuous kinetics, the pathology gradient, the discrimination
#             (AUC) panel and the nested likelihood-ratio panel.
#
# Inputs    : ./kinetics.csv                  (from 01_build_nodule_kinetics_table.py)
#             ./joint_auc_vs_threegroups.csv  (from 02_joint_auc_nested_likelihood_ratio.py)
#
# Outputs   : ./figure1_panels/*.png and *.pdf (one pair per panel);
#             ./figure1_panels/Fig1f_statistics.csv   (written by 04);
#             ./figure1_panels/figure1_panels_summary.md
#
# Usage     : PROJECT_ROOT=/path/to/data python3 03_figure1_panels.py
#             Run 01_build_nodule_kinetics_table.py and
#             02_joint_auc_nested_likelihood_ratio.py first: this script reads the
#             CSV files they leave in this folder.
# -----------------------------------------------------------------------------

import os, math, csv, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Patch
from scipy import stats
from sklearn.metrics import roc_auc_score
from collections import defaultdict, Counter
import statsmodels.api as sm

BASE = os.path.dirname(os.path.abspath(__file__)) + "/"
OUT = BASE + "figure1_panels/"
os.makedirs(OUT, exist_ok=True)

GROUPS = ['Stable', 'GGO-growth', 'SOLID-growth']
# Colours: ggthemes Tableau 10, first three.
COLORS = {'Stable': '#4E79A7', 'GGO-growth': '#F28E2B', 'SOLID-growth': '#E15759'}
GCODE = {'Stable': 0, 'GGO-growth': 1, 'SOLID-growth': 2}
GRADE_LAB = ['Benign', 'AAH', 'AIS', 'MIA', 'IAC']
GRADE_COL = ['#D9D9D9', '#BDD7E7', '#8FBEE0', '#5B9BD5', '#1F4E79']
GREY = '#8EA9C1'
# Binary nodule-type colours: pastel in the Sankey, solid in the composition bars.
BIN_COL_SANKY = {'pGGN': '#FFD7C5', 'PSN': '#C496B8'}
BIN_COL_BAR = {'pGGN': '#F28E2B', 'PSN': '#E15759'}
SANKY_GRP = {'Stable': '#FFF4D1', 'GGO-growth': '#F1DCD4', 'SOLID-growth': '#DCE8E2'}

rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'],
                 'font.size': 8.5, 'axes.linewidth': 0.7, 'xtick.direction': 'out',
                 'ytick.direction': 'out', 'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
                 'figure.dpi': 300, 'savefig.bbox': 'tight'})

LOG = []
def log(s=''):
    print(s); LOG.append(str(s))

def save(fig, name):
    fig.savefig(OUT + name + ".png", dpi=300)
    fig.savefig(OUT + name + ".pdf")
    plt.close(fig)
    log(f"  saved {name}.png / .pdf")

def panel_letter(ax, L, dx=-0.16, dy=1.06):
    ax.text(dx, dy, L, transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom')

def caption_block(fig, text, fs=6.3, pad_in=0.07, color='#333333'):
    """Press the whole caption block against the bottom of the figure, so that
    savefig(bbox='tight') pulls it into the canvas.

    Below the panels there are already several stacked layers of tick numbers,
    axis titles and legends; hard-coding a y coordinate would collide with them.
    So measure the current tight-bbox bottom edge first, then place the text just
    under it -- however many layers lie below, the caption lands on the last line.

    The wrap width is measured, never guessed: panels run from 3.8 in to 7.3 in
    wide, so a character-count estimate is unreliable. Lay out one version, draw
    once to measure the real width, and if it overflows narrow it proportionally
    (at most three attempts).
    """
    fig.canvas.draw()
    bb = fig.get_tightbbox(fig.canvas.get_renderer())        # unit: inch
    # Characters per line: bbox width (in) -> pt, divided by glyph width (pt).
    # Arial averages about 0.48 em; estimate deliberately narrow so the first pass
    # runs wide and the loop below can narrow it -- an over-wide estimate never shrinks.
    ncol = max(28, int(bb.width * 72.0 / (fs * 0.48)))
    fx, fy = fig.get_figwidth(), fig.get_figheight()
    t = None
    for _ in range(4):
        if t is not None:
            t.remove()
        t = fig.text(bb.x0 / fx, (bb.y0 - pad_in) / fy,
                     "\n".join(textwrap.wrap(text, ncol)), ha='left', va='top',
                     fontsize=fs, color=color, linespacing=1.42)
        fig.canvas.draw()
        w = t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
        if abs(w - bb.width) / bb.width < 0.03:
            break
        ncol = max(28, int(ncol * bb.width / w * 0.995))
    # A thin separator rule breaking the caption off from the figure above
    fig.add_artist(plt.Line2D([bb.x0 / fx, (bb.x0 + bb.width) / fx],
                              [(bb.y0 - pad_in * 0.45) / fy] * 2,
                              color='#D0D0D0', lw=0.6, transform=fig.transFigure))
    return t

# ---- In-figure captions: printed at the bottom so each panel explains itself ----
CAP_G = ("Histological grade of the 261 resected nodules across the three progression patterns "
         "(Stable 106, GGO-growth 100, SOLID-growth 55). Numbers inside the bars are nodule "
         "counts. The invasive adenocarcinoma (IAC) share rises monotonically from 25% to 54% "
         "to 91% (Jonckheere–Terpstra z = +7.10, P = 1.3e-12).")
CAP_J = ("Information the three progression patterns add beyond each single variable. For every "
         "metric, a logistic model containing that metric was refitted with the three-pattern "
         "classification added; bars give the likelihood-ratio Δχ² (2 df). Bar colour marks the "
         "denominator available for that metric. All nine metrics sit far past the P = 0.05 "
         "threshold (dashed), the smallest being Δ diameter at Δχ² = 39.9 (P = 2.1e-09).")

# ------------------------------------------------------------ statistics helpers
def num(x):
    try:
        v = float(x)
        return None if (isinstance(x, str) and x.strip() == "") else v
    except (TypeError, ValueError):
        return None

def clean_group(g):
    s = str(g).strip()
    return 'GGO-growth' if s == 'GK253GO-growth' else s

def diam(V):
    return (6.0 * V / math.pi) ** (1.0 / 3.0)

def patho_int(p):
    s = str(p).strip()
    return int(float(s)) if s in ('0', '1', '2', '3', '4') else None

def jonckheere(group_vals, ordered):
    n = [len(group_vals[g]) for g in ordered]
    J = 0.0
    for i in range(len(ordered) - 1):
        for j in range(i + 1, len(ordered)):
            for a in group_vals[ordered[i]]:
                for b in group_vals[ordered[j]]:
                    if a < b:   J += 1.0
                    elif a == b: J += 0.5
    N = sum(n)
    E = (N * N - sum(ni * ni for ni in n)) / 4.0
    Var = (N * N * (2 * N + 3) - sum(ni * ni * (2 * ni + 3) for ni in n)) / 72.0
    z = (J - E) / math.sqrt(Var) if Var > 0 else 0.0
    return z, 2 * (1 - stats.norm.cdf(abs(z)))

def jt_from_labels(labels, patho):
    gv = defaultdict(list)
    for lab, pv in zip(labels, patho):
        gv[lab].append(pv)
    ordered = sorted(gv.keys())
    if len(ordered) < 2 or any(len(gv[o]) == 0 for o in ordered):
        return None
    z, p = jonckheere(gv, ordered)
    grad = [100.0 * sum(1 for x in gv[o] if x == 4) / len(gv[o]) for o in ordered]
    return z, p, grad, ordered

def bh(ps):
    ps = np.asarray(ps, float); k = len(ps)
    order = np.argsort(ps); adj = np.empty(k)
    adj[order] = np.minimum.accumulate((ps[order] * k / np.arange(1, k + 1))[::-1])[::-1]
    return np.clip(adj, 0, 1)

def star(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'

def group_dummies(g):
    return np.column_stack([(g == 1).astype(float), (g == 2).astype(float)])

def lrt(iac, X0, X1):
    m0 = sm.Logit(iac, sm.add_constant(X0)).fit(disp=0)
    m1 = sm.Logit(iac, sm.add_constant(X1)).fit(disp=0)
    df = X1.shape[1] - X0.shape[1]
    chisq = 2 * (m1.llf - m0.llf)
    return chisq, stats.chi2.sf(chisq, df), df

# --------------------------------------------------------- volume conventions
# Two size conventions are used side by side:
#
#   3D convention (vol_pair): delta diameter / relative volume change / annualized
#       relative volume change, from the 3D quantification columns (211/261 nodules,
#       against 147/261 for the manual 2D columns).
#
#   2D convention (bvol/fvol): the Fleischner >= 2 mm criterion, kept on the 2D
#       measurement it is defined on rather than converted to 3D. The comparator is
#       deliberately left on the guideline's own convention.
def vol_pair(r):
    """(first, last) 3D total volume; (None, None) if either is missing or non-positive."""
    b, f = num(r['V0_3d']), num(r['V1_3d'])
    return (b, f) if (b and f and b > 0 and f > 0) else (None, None)

def fleischner_size(r):
    b, f = num(r['bvol']), num(r['fvol'])
    if b and f and b > 0 and f > 0:
        return diam(f) - diam(b) >= 2.0
    return False

def new_solid(r):
    return str(r['btype']) in ('0', '1') and str(r['ftype']) == '2'

# -------------------------------------------------------------------- load
K = list(csv.DictReader(open(BASE + "kinetics.csv", encoding="utf-8-sig")))
for r in K:
    r['group'] = clean_group(r['group'])
K = [r for r in K if r['group'] in GROUPS]
CUR_MODE = K[0].get('MODE', '?') if K else '?'
N_TOT = Counter(r['group'] for r in K)
log(f"kinetics.csv: {len(K)} rows, MODE={CUR_MODE}")
log(f"group totals: " + "  ".join(f"{g}={N_TOT[g]}" for g in GROUPS))

precs = [r for r in K if patho_int(r.get('patho')) is not None]
iac_of = lambda p: 1 if p == 4 else 0


# ======================================================================
# Fig 1c: Sankey -- first-to-last pGGN/PSN transition, one column per group (binary)
# ======================================================================
log("")
log("=" * 78)
log("[c] Sankey: first-to-last transition (binary)")
log("=" * 78)
BIN_LAB = ['pGGN', 'PSN']

def binf(v):
    return 0 if str(v) == '0' else 1

trans = {}
for g in GROUPS:
    d = [r for r in K if r['group'] == g]
    c = np.zeros((2, 2), int)
    for r in d:
        c[binf(r['btype']), binf(r['ftype'])] += 1
    trans[g] = c
    log(f"  {g:14s} n={len(d):4d}  " +
        "  ".join(f"{BIN_LAB[i]}→{BIN_LAB[j]}:{c[i,j]}" for i in range(2) for j in range(2)))

fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.9))
for gi, g in enumerate(GROUPS):
    ax = axes[gi]
    c = trans[g]; n = c.sum()
    ax.set_xlim(-0.15, 1.15); ax.set_ylim(-0.06, 1.10); ax.axis('off')
    # node heights in proportion
    gap, H = 0.045, 0.80
    hn0 = c.sum(axis=1) / n * H
    hn1 = c.sum(axis=0) / n * H
    y0 = np.cumsum([1.00] + list(-hn0))[:-1] - gap
    y1 = np.cumsum([1.00] + list(-hn1))[:-1] - gap
    # nodes
    for i in range(2):
        ax.add_patch(Rectangle((0.0, y0[i] - hn0[i]), 0.07, hn0[i],
                               facecolor=BIN_COL_SANKY[BIN_LAB[i]], edgecolor='none'))
        ax.text(-0.03, y0[i] - hn0[i] / 2, f"{BIN_LAB[i]}\n{int(c.sum(axis=1)[i])}",
                ha='right', va='center', fontsize=6.6, linespacing=1.3)
    for j in range(2):
        ax.add_patch(Rectangle((0.93, y1[j] - hn1[j]), 0.07, hn1[j],
                               facecolor=BIN_COL_SANKY[BIN_LAB[j]], edgecolor='none'))
        ax.text(1.03, y1[j] - hn1[j] / 2, f"{BIN_LAB[j]}\n{int(c.sum(axis=0)[j])}",
                ha='left', va='center', fontsize=6.6, linespacing=1.3)
    # ribbons
    off0 = y0.copy(); off1 = y1.copy()
    for i in range(2):
        for j in range(2):
            v = c[i, j]
            if v == 0:
                continue
            h = v / n * H
            ya, yb = off0[i] - h, off1[j] - h
            t = np.linspace(0, 1, 60)
            s = 3 * t ** 2 - 2 * t ** 3
            xs = 0.07 + (0.93 - 0.07) * t
            top = ya + (yb - ya) * s
            ax.fill_between(xs, top, top + h, color=BIN_COL_SANKY[BIN_LAB[j]], alpha=0.85, lw=0)
            off0[i] -= h; off1[j] -= h
    ax.text(0.5, 1.055, f"{g}  (n = {n})", ha='center', va='bottom', fontsize=7.6,
            color=COLORS[g], fontweight='bold')
    ax.plot([0.035, 0.965], [-0.035, -0.035], color='#BBBBBB', lw=0.6)
    ax.text(0.035, -0.055, "First CT", ha='left', va='top', fontsize=6.6, color='#777777')
    ax.text(0.965, -0.055, "Last CT", ha='right', va='top', fontsize=6.6, color='#777777')
axes[0].text(-0.03, 1.055, 'c', transform=axes[0].transAxes, fontsize=12, fontweight='bold',
             va='bottom', ha='right')
save(fig, "Fig1c_sankey")

# ======================================================================
# Fig 1d: nodule-type composition at the first CT (binary)
# ======================================================================
log("")
log("=" * 78)
log("[d] nodule-type composition at the first CT")
log("=" * 78)
fig, ax = plt.subplots(figsize=(3.5, 3.0))
bottom = np.zeros(3)
comp = {}
for j, lab in enumerate(BIN_LAB):
    vals = []
    for g in GROUPS:
        d = [r for r in K if r['group'] == g]
        vals.append(sum(1 for r in d if binf(r['btype']) == j) / len(d) * 100)
    comp[lab] = np.array(vals)
    ax.bar(range(3), vals, bottom=bottom, color=BIN_COL_BAR[lab], label=lab,
           edgecolor='white', lw=0.6, width=0.62)
    for k, v in enumerate(vals):
        ax.text(k, bottom[k] + v / 2, f"{v:.1f}%", ha='center', va='center',
                fontsize=7.2, color='white' if lab == 'PSN' else '#333333')
    bottom += np.array(vals)
    log(f"  {lab}: " + "  ".join(f"{g}={v:.1f}%" for g, v in zip(GROUPS, vals)))
ax.set_xticks(range(3))
ax.set_xticklabels([f"{g}\n(n = {N_TOT[g]})" for g in GROUPS], fontsize=7.2, linespacing=1.4)
ax.set_ylabel('Nodules at first CT (%)', fontsize=8.5)
ax.set_ylim(0, 100); ax.set_xlim(-0.6, 2.6)
ax.legend(frameon=False, fontsize=7.4, loc='upper center',
          bbox_to_anchor=(0.5, -0.26), ncol=2, handlelength=1.1, handletextpad=0.5)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
panel_letter(ax, 'd', dx=-0.30, dy=1.04)
save(fig, "Fig1d_baseline_composition")


# ======================================================================
# Fig 1f: three-group continuous kinetics (VDT / MDT total / MDT solid / ann. CTR / follow-up)
#   The five-panel version with the CTR display window lives in
#   04_figure1_continuous_kinetics.py so that the two copies cannot drift; this
#   script simply calls it.
# ======================================================================
log("")
log("=" * 78)
log("[f] three-group continuous kinetics")
log("=" * 78)
import subprocess, sys
subprocess.run([sys.executable, BASE + "04_figure1_continuous_kinetics.py"], check=True)
# The statistics are computed on the other side; read them back so this log carries
# the [f] section as well.
_f16 = OUT + "Fig1f_statistics.csv"
if os.path.exists(_f16):
    for _row in csv.DictReader(open(_f16, encoding="utf-8-sig")):
        log(f"  {_row['metric']:14s} KW P={float(_row['kw_p']):.2e}  "
            f"n={_row['n_analyzable']}")
    log("  saved Fig1f_continuous_kinetics.png / .pdf")

# ======================================================================
# Fig 1g: pathology grade 0-4 across the three groups
# ======================================================================
log("")
log("=" * 78)
log("[g] pathology invasiveness spectrum")
log("=" * 78)
cross = {g: Counter() for g in GROUPS}
for r in precs:
    cross[r['group']][patho_int(r['patho'])] += 1
grp_n = {g: sum(cross[g].values()) for g in GROUPS}
patho_groups = {g: [patho_int(r['patho']) for r in precs if r['group'] == g] for g in GROUPS}
z_our, p_our = jonckheere(patho_groups, GROUPS)
log(f"  evaluable pathology n={len(precs)}  " + "  ".join(f"{g}={grp_n[g]}" for g in GROUPS))
for g in GROUPS:
    log(f"  {g:14s} " + "  ".join(f"{GRADE_LAB[j]}:{cross[g][j]}" for j in range(5)) +
        f"   IAC={cross[g][4]/grp_n[g]*100:.1f}%")
log(f"  JT z={z_our:+.2f}  P={p_our:.2e}")

fig, ax = plt.subplots(figsize=(3.6, 3.6))
bottom = np.zeros(3)
for j in range(5):
    vals = [cross[g][j] / grp_n[g] * 100 for g in GROUPS]
    ax.bar(range(3), vals, bottom=bottom, color=GRADE_COL[j], label=GRADE_LAB[j],
           edgecolor='white', lw=0.5, width=0.62)
    for k, v in enumerate(vals):
        if v >= 6:
            ax.text(k, bottom[k] + v / 2, f"{v:.0f}", ha='center', va='center',
                    fontsize=6.6, color='white' if j >= 3 else '#333333')
    bottom += np.array(vals)
ax.set_xticks(range(3))
ax.set_xticklabels([f"{g}\n(n = {grp_n[g]})" for g in GROUPS], fontsize=7.2, linespacing=1.4)
ax.set_ylabel('Resected nodules (%)', fontsize=8.5)
ax.set_ylim(0, 100); ax.set_xlim(-0.6, 2.6)
# The JT annotation used to sit at (0.97, 0.97) -- inside the dark-blue SOLID-growth
# bar, where grey on navy is unreadable and squashes the in-bar "89" to the edge.
# Moved outside the axes at the top, leaving the inside of the bars to the values.
ax.text(0.99, 1.012, f'JT z = {z_our:+.2f}    P = {p_our:.1e}', transform=ax.transAxes,
        ha='right', va='bottom', fontsize=6.8, color='#333333')
# With ncol=3 matplotlib fills column-first, splitting the five-grade ramp into
# [Benign, AIS, IAC] / [AAH, MIA] -- the reading order scrambles and the gradient
# breaks. One row of five keeps the colour ramp continuous.
ax.legend(frameon=False, fontsize=6.2, loc='upper center', bbox_to_anchor=(0.5, -0.17),
          ncol=5, columnspacing=0.7, handlelength=0.85, handletextpad=0.35)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
panel_letter(ax, 'g', dx=-0.28, dy=1.04)
caption_block(fig, CAP_G)
save(fig, "Fig1g_pathology_gradient")


# ======================================================================
# Fig 1i: discrimination (AUC)
# ======================================================================
log("")
log("=" * 78)
log("[i] discrimination (AUC)")
log("=" * 78)
our_labels_full = np.array([GCODE[r['group']] for r in precs])
our_patho_full = np.array([patho_int(r['patho']) for r in precs])
our_iac_full = np.array([iac_of(patho_int(r['patho'])) for r in precs])
our_auc_full = roc_auc_score(our_iac_full, our_labels_full)
log(f"  Our 3 groups, all pathology (n={len(precs)}): AUC = {our_auc_full:.3f}")

def f_dd(r):
    b, f = vol_pair(r)
    return (diam(f) - diam(b)) if (b and f) else None

def f_rel(r):
    b, f = vol_pair(r)
    return ((f - b) / b) if (b and f) else None

def f_annrel(r):
    rel, fu = f_rel(r), num(r['fu'])
    return (rel / (fu / 365.25)) if (rel is not None and fu and fu > 0) else None

def f_vdt(r):
    v = num(r['VDT'])
    return v if (v is not None and v > 0) else None

def f_mdt(r):
    v = num(r['MDT3d_total'])
    return v if (v is not None and v > 0) else None

def f_ctr(r):
    v = num(r.get('ann_CTR3d'))
    return v if (v is not None and np.isfinite(v)) else None

AUCVARS = [
    ('VDT',                          f_vdt,    -1),
    ('MDT, total mass',              f_mdt,    -1),
    ('Annualized CTR change',        f_ctr,    +1),
    ('Δ diameter (mm)',              f_dd,     +1),
    ('Annualized rel. vol. change',  f_annrel, +1),
    ('Relative volume change',       f_rel,    +1),
]
auc_rows = []
for name, fn, direction in AUCVARS:
    vals, iac, grank = [], [], []
    for r in precs:
        v = fn(r)
        if v is None or not np.isfinite(v):
            continue
        vals.append(v)
        iac.append(iac_of(patho_int(r['patho'])))
        grank.append(GCODE[r['group']])
    vals, iac, grank = np.array(vals), np.array(iac), np.array(grank)
    if len(vals) < 20 or len(set(iac)) < 2:
        continue
    auc_c = roc_auc_score(iac, vals * direction)
    auc_g = roc_auc_score(iac, grank)
    auc_rows.append((name, auc_c, auc_g, len(vals)))
    log(f"  {name:30s} n={len(vals):4d}  AUC(continuous)={auc_c:.3f}  AUC(3 groups, same n)={auc_g:.3f}")

# ----------------------------------------------------------------------
# Layout: bands by "denominator" + an in-row dumbbell
#   The six metrics each have their own n; the old side-by-side bars made the reader
#   count n to know which rows are comparable. The bands below group rows that share
#   an n and name the denominator directly.
#   In-row dumbbell: open grey dot = the single metric's own AUC, filled red dot =
#   the three groups on *the same nodules*; the connecting line is the extra
#   discrimination the grouping adds over that continuous metric.
#   The six red dots fall in 0.732-0.736, almost a vertical line -- the answer does
#   not move with the denominator, the direct reply to "the n differ, so the AUCs are
#   not comparable", without turning to the supplement.
#
#   Row tuple = (name, n, left AUC, right AUC, left marker). Left marker 'o' = an
#   unfitted rank statistic, 'D' = a fitted multi-variable logistic (which requires
#   the CV AUC, see the variant below). The "Our 3 groups" row has None on the right
#   and takes the solid-bar branch, so it does not use the marker slot.
# ----------------------------------------------------------------------
def draw_fig1i(bands, cap, fname, mlegend):
    """Draw one discrimination panel. bands sets how many bands there are and how
    many rows each holds.

    The height adapts to the band count: the row pitch is in fixed axis units, so more
    bands need a taller canvas -- otherwise each row is squashed and the two-line
    labels collide.
    """
    span = sum(0.76 + len(items) + 0.36 for _, items in bands)
    fig, ax = plt.subplots(figsize=(6.9, span * (3.9 / 10.36)))
    XMIN, XMAX = 0.472, 0.792

    # ---- vertical layout: band header 0.76, row pitch 1.00, inter-band gap 0.36 ----
    placed, spans, y = [], [], 0.0
    for bi, (bname, items) in enumerate(bands):
        y_hdr = y; y -= 0.76
        for nm, n, c, g, mk in items:
            placed.append((bi, nm, n, c, g, mk, y)); y -= 1.0
        spans.append((bi, bname, y + 1.0, y_hdr))
        y -= 0.36

    # ---- band shading: alternate light grey, drawing "these rows share a denominator" ----
    for bi, bname, ylo, yhi in spans:
        if bi % 2 == 0:
            ax.axhspan(ylo - 0.54, yhi + 0.30, color='#F4F6F8', zorder=0)

    # ---- band name: inside the canvas at the top-left, small caps, so it does not fight
    # the off-axis labels. White bbox: the vertical chance line would pass through it;
    # a backing layer pushes the line behind the text ----
    for bi, bname, ylo, yhi in spans:
        ax.text(0.008, yhi, bname.upper(), ha='left', va='center', fontsize=6.2,
                color='#9C9C9C', fontweight='bold', zorder=5,
                bbox=dict(facecolor='white', edgecolor='none', pad=1.0, alpha=0.92),
                transform=ax.get_yaxis_transform())

    # ---- row labels: name + n on two lines, n alone so rows align for comparison ----
    ax.set_yticks([p[6] for p in placed])
    ax.set_yticklabels([f"{p[1]}\nn = {p[2]}" for p in placed], fontsize=6.9, linespacing=1.45)
    for _t, p in zip(ax.get_yticklabels(), placed):
        _top = p[4] is None
        _t.set_color(COLORS['SOLID-growth'] if _top else '#333333')
        _t.set_fontweight('bold' if _top else 'normal')
    ax.tick_params(axis='y', length=0)

    # ---- the top row is a solid bar, the rest are dumbbells ----
    for bi, nm, n, c, g, mk, yy in placed:
        if g is None:
            ax.barh(yy, our_auc_full - XMIN, left=XMIN, height=0.50,
                    color=COLORS['SOLID-growth'], edgecolor='none', zorder=3)
            ax.text(our_auc_full - 0.007, yy, f"{our_auc_full:.3f}", va='center', ha='right',
                    fontsize=8.2, fontweight='bold', color='white', zorder=4)
        else:
            ax.plot([c, g], [yy, yy], color='#C8C8C8', lw=1.7, solid_capstyle='round', zorder=1)
            if mk == 'D':
                # fitted model: a diamond with a dark-grey edge, distinct in both shape
                # and colour from the unfitted single-metric circle
                ax.scatter([c], [yy], s=30, marker='D', facecolor='white',
                           edgecolor='#555555', linewidths=1.4, zorder=3)
            else:
                ax.scatter([c], [yy], s=33, facecolor='white', edgecolor=GREY,
                           linewidths=1.5, zorder=3)
            ax.scatter([g], [yy], s=40, color=COLORS['SOLID-growth'], edgecolor='white',
                       linewidths=0.8, zorder=3)
            # white bbox: the grey dot's value can land on the chance line; a backing
            # layer pushes the line behind the text
            ax.text(c - 0.007, yy, f"{c:.3f}", va='center', ha='right', fontsize=6.7,
                    color='#7A7A7A', zorder=4,
                    bbox=dict(facecolor='white', edgecolor='none', pad=0.7, alpha=0.92))
            ax.text(g + 0.007, yy, f"{g:.3f}", va='center', ha='left', fontsize=6.9,
                    color=COLORS['SOLID-growth'], fontweight='bold', zorder=4,
                    bbox=dict(facecolor='white', edgecolor='none', pad=0.7, alpha=0.92))

    ax.axvline(0.5, color='#333333', lw=0.8, ls='--', zorder=1)
    ax.set_xlim(XMIN, XMAX)
    _ytop = max(p[6] for p in placed) + 0.60
    _ybot = min(p[6] for p in placed) - 0.60
    ax.set_ylim(_ybot, _ytop)
    ax.text(0.5, _ybot + 0.04, 'chance', ha='center', va='bottom', fontsize=6.4, color='#777777',
            bbox=dict(facecolor='white', edgecolor='none', pad=0.6, alpha=0.9), zorder=5)
    # What the two dumbbell ends mean goes on the line under the axis title -- the red
    # dots forming one vertical line is this panel's main evidence, and a reader must
    # see at once that they are "the three groups on the same nodules".
    ax.set_xlabel('AUC  (invasive adenocarcinoma vs. rest)\n' + mlegend,
                  fontsize=8.2, linespacing=2.0)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    panel_letter(ax, 'i', dx=-0.075, dy=1.01)
    caption_block(fig, cap)
    save(fig, fname)


BANDS = [
    ('All pathology', [('Our 3 groups', len(precs), our_auc_full, None, 'o')]),
    # Bands follow each metric's *domain*, not an n threshold -- n moves with the
    # convention and a threshold would misplace rows. Rows within a band may have
    # different n (annualized CTR 239 and delta diameter 211 both cover "all nodules",
    # but CTR needs a 3D quantification while diameter needs two 2D volumes); each row
    # carries its own n label.
    ('Growing nodules only',
     [(nm, n, c, g, 'o') for nm, c, g, n in auc_rows if nm.startswith(('VDT', 'MDT'))]),
    ('All nodules with serial volumetric quantification',
     [(nm, n, c, g, 'o') for nm, c, g, n in auc_rows if not nm.startswith(('VDT', 'MDT'))]),
]


# ----------------------------------------------------------------------
# Discrimination band, multi-variable leg
#   A single metric is an unfitted rank statistic while a combination model fits
#   coefficients, so both ends of the band use the repeated 5-fold CV AUC on the
#   same folds (same n, same folds, standardisation inside the training fold only).
#   The values are read from joint_auc_vs_threegroups.csv.
# ----------------------------------------------------------------------
JOINT_LAB = {
    'VDT + MDT + ann_CTR': 'VDT + MDT + CTR',
    'MDT + ann_CTR':       'MDT + CTR',
    'VDT + ann_CTR':       'VDT + CTR',
    'VDT + MDT':           'VDT + MDT',
}
joint_items = []
with open(BASE + "joint_auc_vs_threegroups.csv", encoding="utf-8-sig") as f:
    for d in csv.DictReader(f):
        if int(d['k']) < 2 or d['label'] not in JOINT_LAB:
            continue
        joint_items.append((JOINT_LAB[d['label']], int(d['n']),
                            float(d['cv']), float(d['cv_g']), 'D'))
joint_items.sort(key=lambda r: -r[2])          # strongest combination on top
for _it in joint_items:
    log(f"  [variant] {_it[0]:18s} n={_it[1]:4d}  CV AUC={_it[2]:.3f}  "
        f"three groups, same folds CV={_it[3]:.3f}")

CAP_I_JOINT = ("Open circles are the AUC of each continuous metric on its own; open diamonds are "
               "multi-variable logistic models. Every model is compared with the three progression "
               "patterns evaluated on exactly the same nodules (filled circles), so each pair shares "
               "one denominator. Rows are grouped by where each metric is defined. Single metrics are "
               "unfitted rank statistics; multi-variable models are fitted and are therefore reported "
               "as repeated 5-fold cross-validated AUC, with the three-pattern comparator "
               "cross-validated on the same folds. CTR denotes the annualized CTR change. Across all "
               "denominators the three-pattern AUC stays between 0.736 and 0.746; no single metric "
               "exceeds 0.624 and no multi-variable combination exceeds 0.676.")

BANDS_J = BANDS + [('All combinations of the same variables', joint_items)]
# The legend uses the three marker glyphs Arial actually has: U+25CB / U+25CA / U+25CF.
# Do not use U+25C6 -- Arial lacks that code point and it renders as a tofu box.
draw_fig1i(BANDS_J, CAP_I_JOINT, "Fig1i_discrimination_joint",
           '○ single continuous metric   ◊ fitted multi-variable logistic (5-fold CV)   '
           '● our 3 groups, same subset')

# ======================================================================
# Fig 1j: nested LR -- information the three groups add on top of each single variable
# ======================================================================
log("")
log("=" * 78)
log("[j] nested LR independent information")
log("=" * 78)
iac_full_a = np.array([iac_of(patho_int(r['patho'])) for r in precs])
gr_full = np.array([GCODE[r['group']] for r in precs])
G_full = group_dummies(gr_full)
lr_rows = []
for label, fn in [('New solid component (densification)', new_solid),
                  ('Fleischner ≥2 mm (size)', fleischner_size),
                  ('Fleischner growth (size or densification)',
                   lambda r: fleischner_size(r) or new_solid(r))]:
    X = np.array([int(fn(r)) for r in precs])[:, None]
    chisq, pval, _ = lrt(iac_full_a, X, np.hstack([X, G_full]))
    lr_rows.append((label, chisq, pval, len(precs), True))

for name, fn, direction in AUCVARS:
    vals, iac, grank = [], [], []
    for r in precs:
        v = fn(r)
        if v is None or not np.isfinite(v):
            continue
        vals.append(v); iac.append(iac_of(patho_int(r['patho']))); grank.append(GCODE[r['group']])
    vals, iac, grank = np.array(vals), np.array(iac), np.array(grank)
    if len(vals) < 20:
        continue
    X = (vals * direction)[:, None]
    chisq, pval, _ = lrt(iac, X, np.hstack([X, group_dummies(grank)]))
    lr_rows.append((name, chisq, pval, len(vals), False))

lr_rows.sort(key=lambda r: r[1], reverse=True)   # strongest on top
for label, chisq, pval, n, is_full in lr_rows:
    log(f"  {label:46s} n={n:4d}  Δχ²={chisq:5.1f}  P={pval:.1e}")

# ----------------------------------------------------------------------
# Layout: colour by denominator
#   Bars are ordered by decreasing delta chi-square (strongest on top) and coloured by
#   the analysis set the row used:
#     red    = all pathology, n = 261 (guideline criteria, no volumetric measurement)
#     blue   = all nodules with 3D quantification, n = 239 / 211
#     orange = growing nodules only, n = 206 / 177
#   The colour is assigned from each metric's domain, not from a threshold on n.
# ----------------------------------------------------------------------
C_FULL, C_VOL, C_GROW = '#E15759', '#4E79A7', '#F28E2B'


def _ncol(name, n):
    if name == 'Our 3 groups' or n >= 261:     # all pathology: guideline criteria, no volumetry
        return C_FULL
    if name.startswith(('VDT', 'MDT')):        # defined only for growing nodules
        return C_GROW
    return C_VOL

fig, ax = plt.subplots(figsize=(6.4, 3.6))
ys = np.arange(len(lr_rows))[::-1]
cs = [r[1] for r in lr_rows]
ax.barh(ys, cs, color=[_ncol(r[0], r[3]) for r in lr_rows], edgecolor='white',
        height=0.68, zorder=3)
for y, (label, chisq, pval, n, is_full) in zip(ys, lr_rows):
    ax.text(chisq + max(cs) * 0.014, y + 0.17, f"{chisq:.1f}", va='center', fontsize=7.8,
            color=_ncol(label, n), fontweight='bold', zorder=4)
    ax.text(chisq + max(cs) * 0.014, y - 0.23, f"P = {pval:.1e}", va='center', fontsize=6.3,
            color='#999999', zorder=4)
ax.set_yticks(ys)
ax.set_yticklabels([f"{r[0]}\nn = {r[3]}" for r in lr_rows], fontsize=7.0, linespacing=1.40)
ax.tick_params(axis='y', length=0)
ax.set_xlabel('Δχ²  (likelihood-ratio, 2 df) — information the 3 groups add\n'
              'beyond the variable, in a model that already contains it',
              fontsize=7.8, linespacing=1.5)
ax.set_xlim(0, max(cs) * 1.26)
# The threshold label used to sit vertically inside the bars, grey on colour and
# unreadable; raise it above every bar, horizontally.
ax.set_ylim(-0.78, len(lr_rows) + 0.32)
_chi = stats.chi2.ppf(0.95, 2)
ax.axvline(_chi, color='#333333', lw=0.8, ls='--', zorder=2)
ax.text(_chi + max(cs) * 0.005, len(lr_rows) - 0.26, 'P = 0.05', fontsize=6.3,
        color='#777777', ha='left', va='center')
# denominator legend: each colour says how many nodules the row used
ax.legend(handles=[Patch(facecolor=C_FULL, edgecolor='white',
                         label='all pathology  (n = 261)'),
                   Patch(facecolor=C_VOL, edgecolor='white',
                         label='serial volumes  (n = 239 / 211)'),
                   Patch(facecolor=C_GROW, edgecolor='white',
                         label='growing only  (n = 206 / 177)')],
          frameon=False, fontsize=6.3, loc='center right',
          bbox_to_anchor=(1.005, 0.33), handlelength=0.9, handletextpad=0.5,
          labelspacing=0.35, borderpad=0.0)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
panel_letter(ax, 'j', dx=-0.44, dy=1.01)
caption_block(fig, CAP_J)
save(fig, "Fig1j_independent_information")

# ------------------------------------------------------------------- summary
log("")
log("=" * 78)
log(f"Done. Panel output directory: {OUT}")
log("=" * 78)
with open(OUT + "figure1_panels_summary.md", "w", encoding="utf-8") as f:
    f.write("# Clinical-kinetics panels -- results\n\n")
    f.write(f"- Data source: `kinetics.csv` ({len(K)} rows, 999 nodules / 845 patients)\n")
    f.write(f"- Primary VDT convention = `{CUR_MODE}`\n")
    f.write("- Pathology 0 benign < 1 AAH < 2 AIS < 3 MIA < 4 IAC; binary nodule type 0 = pGGN, 1-2 = PSN\n\n")
    f.write("```\n" + "\n".join(LOG) + "\n```\n")
log("  summary written to figure1_panels_summary.md")
