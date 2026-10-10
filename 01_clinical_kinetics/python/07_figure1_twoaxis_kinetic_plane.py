#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Two-axis kinetic plane
#
# Purpose   : Draw the two-axis kinetic plane at the final print size: a
#             whole-cohort grey 2-D KDE background with three quantile contours,
#             three-group median/IQR crosses with individual nodules, ordered
#             progression arrows and 1-D density margins with ticks.
#
# Inputs    : $PROJECT_ROOT/omics/supplementary_clinical_data/VDT_MDT_continuous_kinetics_VDT2d_fill3d/kinetics.csv
#
# Outputs   : ./figures/Fig1c_twoaxis.pdf / .png
#
# Usage     : PROJECT_ROOT=/path/to/data python3 07_figure1_twoaxis_kinetic_plane.py
# -----------------------------------------------------------------------------
"""
Two-axis kinetic plane
=============================================================================
Layout: a whole-cohort grey 2-D KDE background with three quantile contours,
three-group median / IQR crosses with individual nodules, ordered progression
arrows, and 1-D density margins with ticks on the top and the right.

Colours: Stable #4E79A7 / GGO-growth #F28E2B / SOLID-growth #E15759.
The display window, the clipping definition and the VDT = 400 d reference line
are the ones used throughout the clinical-kinetics analyses.

Sizing
  figsize is the final print size; bbox_inches='tight' is **not** used (it would
  change the canvas size and the font size could no longer be derived). Place it
  at 1:1 when assembling; what you see is the font size you get.
  Font floor 6.6 pt.

Data source
  VDT_MDT_continuous_kinetics_VDT2d_fill3d/kinetics.csv   (999 nodules x 45 columns)
  Columns are read by header name, never by position, to avoid column mix-ups
=============================================================================
"""
# -----------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", "/path/to/GGN_progression_data")
PROJ = PROJECT_ROOT.rstrip("/") + "/"
VDT = os.path.join(PROJ, "omics/supplementary_clinical_data/VDT_MDT_continuous_kinetics_VDT2d_fill3d")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)

# ---- Definitions -----------------------------------------------------------
GROUPS = ["Stable", "GGO-growth", "SOLID-growth"]
COLORS = {"Stable": "#4E79A7", "GGO-growth": "#F28E2B", "SOLID-growth": "#E15759"}
CLIP_SG, CLIP_SOLID = 0.2, 0.3          # display-window clipping: Stable/GGO |dCTR| <= 0.2, SOLID 0-0.3
VREF = 400.0                            # VDT = 400 d (fast-growth landmark)
GRX = lambda vdt: np.log10(365.25 / vdt)   # growth-rate axis = doublings per year
XL, YL = (-2.60, 1.05), (-CLIP_SG - 0.02, CLIP_SOLID + 0.02)

# ---- Font sizes (pt, print size) -------------------------------------------
FS_TICK, FS_AXIS, FS_GRP, FS_NOTE = 7.0, 7.4, 7.0, 6.6

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,     # fonts embedded as TrueType, a submission requirement
    "axes.linewidth": 0.6,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "figure.dpi": 300,
})

# ============================================================================
# 1. Load the data
# ============================================================================
df = pd.read_csv(os.path.join(VDT, "kinetics.csv"))
df["group"] = df["group"].replace({"GK253GO-growth": "GGO-growth"})
df = df[df["group"].isin(GROUPS)].copy()

pts = {}
for g in GROUPS:
    d = df[df["group"] == g].dropna(subset=["VDT", "ann_CTR3d"])
    d = d[d["VDT"] > 0]
    if g == "SOLID-growth":
        d = d[(d["ann_CTR3d"] >= 0) & (d["ann_CTR3d"] <= CLIP_SOLID)]
    else:
        d = d[d["ann_CTR3d"].abs() <= CLIP_SG]
    pts[g] = (GRX(d["VDT"].to_numpy()), d["ann_CTR3d"].to_numpy())

for g in GROUPS:
    x, y = pts[g]
    print(f"{g:14s} n(display)={len(x):4d}   median VDT={365.25/10**np.median(x):6.0f} d"
          f"   median annCTR={np.median(y):+.4f}/yr")

# ============================================================================
# 2. Canvas: figsize is the print size; the four margins are set precisely with subplots_adjust, no tight cropping
# ============================================================================
fig_w, fig_h = 3.72, 3.08
fig = plt.figure(figsize=(fig_w, fig_h))
gs = fig.add_gridspec(2, 2, width_ratios=[5, 1.15], height_ratios=[1.25, 5],
                      wspace=0.06, hspace=0.06)
ax_top = fig.add_subplot(gs[0, 0])
ax = fig.add_subplot(gs[1, 0])
ax_right = fig.add_subplot(gs[1, 1])

gx = np.linspace(*XL, 220)
gy = np.linspace(*YL, 220)
GX, GY = np.meshgrid(gx, gy)
POS = np.vstack([GX.ravel(), GY.ravel()])

# ============================================================================
# 3. Main plane: whole-cohort grey 2-D kernel density (only three quantile contours, clean)
# ============================================================================
pool_x = np.concatenate([pts[g][0] for g in GROUPS])
pool_y = np.concatenate([pts[g][1] for g in GROUPS])
Z = stats.gaussian_kde(np.vstack([pool_x, pool_y]))(POS).reshape(GX.shape)

zsa = np.sort(Z.ravel())[::-1]
cuma = np.cumsum(zsa) / zsa.sum()
LV = sorted(zsa[np.searchsorted(cuma, q)] for q in (0.25, 0.50, 0.75))

ax.contourf(GX, GY, Z, levels=[LV[0], Z.max() * 1.01], colors=["#BDBDBD"],
            alpha=0.55, zorder=1)
ax.contour(GX, GY, Z, levels=LV, colors=["#7A7A7A"], linewidths=[0.5, 0.8, 0.5],
           linestyles=[":", "-", ":"], zorder=1)

# ---- Three groups: median +/- IQR crosses + centroid dots ------------------
MED = {}
for g in GROUPS:
    x, y = pts[g]
    mx, my = float(np.median(x)), float(np.median(y))
    qx = np.percentile(x, [25, 75])
    qy = np.percentile(y, [25, 75])
    MED[g] = (mx, my)
    ax.plot([qx[0], qx[1]], [my, my], lw=1.0, color=COLORS[g], alpha=0.95,
            solid_capstyle="butt", zorder=5)
    ax.plot([mx, mx], [qy[0], qy[1]], lw=1.0, color=COLORS[g], alpha=0.95,
            solid_capstyle="butt", zorder=5)
    ax.scatter([mx], [my], s=52, marker="o", c=COLORS[g], edgecolors="black",
               linewidths=0.7, zorder=6)

# ---- Ordered progression arrows -------------------------------------------
ax.annotate("", xy=MED["GGO-growth"], xytext=MED["Stable"],
            arrowprops=dict(arrowstyle="-|>", color="#3A3A3A", lw=1.0,
                            shrinkA=7, shrinkB=8, alpha=0.9, mutation_scale=8), zorder=5)
_g, _s = MED["GGO-growth"], MED["SOLID-growth"]
_dx = 0.30 * (XL[1] - XL[0]) / 3.65
ax.annotate("", xy=(_s[0] + _dx, _s[1] - 0.028), xytext=(_g[0] + _dx, _g[1] + 0.010),
            arrowprops=dict(arrowstyle="-|>", color="#3A3A3A", lw=1.0, alpha=0.9,
                            mutation_scale=8, shrinkA=0, shrinkB=0), zorder=5)

LBL_OFF = {"Stable": (-2, 21), "GGO-growth": (-10, -26), "SOLID-growth": (22, 19)}
for g in GROUPS:
    ax.annotate(f"{g}\nn = {len(pts[g][0])}", MED[g], textcoords="offset points",
                xytext=LBL_OFF[g], ha="center", fontsize=FS_GRP, fontweight="bold",
                color=COLORS[g], zorder=7,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))

# ---- Reference lines + axes -----------------------------------------------
ax.axhline(0, ls="--", lw=0.7, c="grey", zorder=4)
ax.axvline(GRX(VREF), ls=":", lw=0.8, c="grey", zorder=4)
ax.set_xlim(*XL)
ax.set_ylim(*YL)

XT = [0.03, 0.1, 0.3, 1, 3]
ax.set_xticks(np.log10(XT))
ax.set_xticklabels([f"{t:g}" for t in XT], fontsize=FS_TICK)

# The x axis is **computed**, not a column of the source table:
# growth rate = log10(365.25 / VDT days). Writing only "doublings per year" leaves
# the reader unable to invert it, so two things are done together:
#   (1) the axis title states the formula; (2) under each tick sits a line with the
# equivalent VDT in days, so that x=1 <-> VDT 365 days is visible at a glance,
# without flipping back to the main text.
_eq = [f"{365.25 / t:,.0f}" for t in XT]          # 12175 / 3652 / 1218 / 365 / 122
_ytxt = YL[0] - (YL[1] - YL[0]) * 0.148           # sits one line below the major tick numbers
for _t, _v in zip(XT, _eq):
    ax.text(np.log10(_t), _ytxt, _v, transform=ax.transData,
            ha="center", va="top", fontsize=FS_NOTE - 0.6, color="#6E6E6E")
ax.text(XL[0] + 0.06, _ytxt, "VDT (d)", transform=ax.transData,
        ha="left", va="top", fontsize=FS_NOTE - 0.6, color="#6E6E6E")

ax.set_xlabel("Growth rate  =  log10(365.25 / VDT)\ndoublings per year,  log scale",
              fontsize=FS_AXIS, labelpad=13, linespacing=1.55)
ax.tick_params(axis="y", labelsize=FS_TICK)
ax.set_ylabel("Annualized CTR change  (/yr)", fontsize=FS_AXIS)

ax.text(GRX(VREF) + 0.04, YL[0] + 0.008, f"VDT = {VREF:.0f} d", rotation=90,
        va="bottom", ha="left", fontsize=FS_NOTE, color="grey")
# Two legends stacked into a two-line key, in the blank at the top-left of the main plane (the middle is filled by the grey density and the group labels)
ax.text(XL[0] + 0.06, YL[1] - 0.012, "grey: pooled density (all nodules)",
        ha="left", va="top", fontsize=FS_NOTE, color="#8C8C8C")
ax.text(XL[0] + 0.06, YL[1] - 0.038, "dashed: no densification (ΔCTR/yr=0)",
        ha="left", va="top", fontsize=FS_NOTE, color="#8C8C8C")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ============================================================================
# 4. Marginal 1-D kernel densities (top offset by group, with ticks; right likewise)
# ============================================================================
xgrid = np.linspace(*XL, 300)
ygrad = np.linspace(*YL, 300)
for gi, g in enumerate(GROUPS):
    x, y = pts[g]
    kx = stats.gaussian_kde(x)(xgrid)
    kx = kx / kx.max() * 0.62
    ax_top.fill_between(xgrid, gi, gi + kx, color=COLORS[g], alpha=0.55, lw=0)
    ax_top.plot(xgrid, gi + kx, color=COLORS[g], lw=0.7)

    ky = stats.gaussian_kde(y)(ygrad)
    ky = ky / ky.max() * 0.62
    ax_right.fill_betweenx(ygrad, gi, gi + ky, color=COLORS[g], alpha=0.55, lw=0)
    ax_right.plot(gi + ky, ygrad, color=COLORS[g], lw=0.7)

ax_top.set_xlim(*XL)
ax_top.set_ylim(-0.05, 3.0)
# The group name sits in the blank at the left end of each ridge - the panel is only 3.6 in wide and cannot hold a long left label
ax_top.set_yticks([])
ax_top.set_xticks(np.log10(XT))
ax_top.tick_params(axis="x", labelbottom=False)
ax_top.axvline(GRX(VREF), ls=":", lw=0.8, c="grey")
for gi, g in enumerate(GROUPS):
    ax_top.text(XL[0] + 0.05, gi + 0.55, g, ha="left", va="center", fontsize=FS_GRP,
                color=COLORS[g], zorder=5,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8))
for s in ["top", "right", "bottom", "left"]:
    ax_top.spines[s].set_visible(False)

ax_right.set_ylim(*YL)
ax_right.set_xlim(-0.05, 3.0)
ax_right.set_xticks([0.5, 1.5, 2.5])
ax_right.set_xticklabels([])          # the colour key is given by the top ridge; here the ticks stay only as a position reference
ax_right.tick_params(axis="x", length=2.2)
ax_right.tick_params(axis="y", labelleft=False)
ax_right.axhline(0, ls="--", lw=0.7, c="grey")
for s in ["top", "left", "right"]:
    ax_right.spines[s].set_visible(False)
ax_right.spines["bottom"].set_linewidth(0.5)

# ============================================================================
# 5. Output (1:1 print size, no tight cropping)
# ============================================================================
# bottom raised from 0.175 to 0.265: three things now sit below the axis - the
# major tick numbers, the VDT-days line and the two-line axis title. 0.175 fits
# only the first two, and the second title line "doublings per year, log scale"
# would be cut off by the bottom edge (only the first line renders). 0.265 of the
# figure height is about 0.82 in, just enough for all three.
fig.subplots_adjust(left=0.128, right=0.975, top=0.975, bottom=0.265)

png = os.path.join(FIG, "Fig1c_twoaxis.png")
fig.savefig(os.path.join(FIG, "Fig1c_twoaxis.pdf"), format="pdf")
fig.savefig(png, dpi=300)
plt.close(fig)

print(f"\ncanvas {fig_w} x {fig_h} in (= print size, no scaling)")
print(f"smallest font {FS_NOTE} pt, axis label {FS_AXIS} pt")
print(f"→ figures/Fig1c_twoaxis.pdf / .png")
