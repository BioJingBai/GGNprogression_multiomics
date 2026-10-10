#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF step 7 - distance from SPP1+ macrophages to exhausted T-cell subsets
#
# Purpose   : Two-panel figure on one shared linear axis: grouped boxplots with
#             jittered points and pairwise brackets, and group-median trajectory
#             lines with a vertical significance marker per subset.
#
# Inputs    : <OUTPUT_DATA>/spatial_results/SPP1Mac_PD1_CD8_NND_median.csv
#             <OUTPUT_DATA>/spatial_results/SPP1Mac_PD1_CD4_NND_median.csv
#             <OUTPUT_DATA>/spatial_results/SPP1Mac_PD1_NND_median.csv
#             all produced by 03_final_spatial.py
#
# Outputs   : <FIGS_DIR>/preview/05_distance_shared_axis.{pdf,png}
#             <FIGS_DIR>/03_Macrophage_Polarization/
#                 Fig3_SPP1Mac_ExhaustedT_Distance_shared_axis.{pdf,png}
#
# Usage     : PROJECT_ROOT=/path/to/data python 07_distance_shared_axis.py
# -----------------------------------------------------------------------------

import os
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import figstyle as fs

PROJECT_ROOT = os.environ.get("PROJECT_ROOT", "/path/to/GGN_progression_data")
BASE = os.path.join(PROJECT_ROOT, "mIF_Validation_Results")
DATA = os.path.join(BASE, "data")
OUT = os.path.join(BASE, "figures", "preview")
OUT_FORMAL = os.path.join(BASE, "figures", "03_Macrophage_Polarization")
NAME_FORMAL = "Fig3_SPP1Mac_ExhaustedT_Distance_shared_axis"
fs.prepare()

GRPS = ["Stable", "Growth", "Solid"]
GLAB, GCOL = fs.GROUP_LABELS, fs.GROUP_COLORS
SUBSETS = [("SPP1Mac_PD1_CD8", "PD-1+CD8+"), ("SPP1Mac_PD1_CD4", "PD-1+CD4+"),
           ("SPP1Mac_PD1", "PD-1+")]
BRACKET_FRAC = [0.885, 0.925, 0.968]   # panel A: place brackets in axis fraction, never stretching the y-axis

VAL = {}
for key, lab in SUBSETS:
    d = pd.read_csv(os.path.join(DATA, "spatial_results", f"{key}_NND_median.csv"))
    VAL[lab] = {g: d.loc[d["group"] == g, "obs_nnd"].dropna().values.astype(float)
                for g in GRPS}


def stats(lab):
    from scipy.stats import kruskal, mannwhitneyu
    d = {g: v for g, v in VAL[lab].items() if len(v) >= 2}
    try:
        kw = kruskal(*d.values())[1]
    except ValueError:
        kw = 1.0
    pairs = []
    for g1, g2 in itertools.combinations(GRPS, 2):
        try:
            p = mannwhitneyu(d[g1], d[g2], alternative="two-sided")[1]
        except (ValueError, KeyError):
            p = np.nan
        pairs.append(((g1, g2), min(p * 3, 1.0) if p == p else np.nan))
    return d, kw, pairs


def fmt(p):
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}"


def nest(lab, pairs):
    """Order the three pairs by group median so the brackets nest (narrow below,
    wide above) instead of crossing."""
    med = {g: float(np.median(VAL[lab][g])) for g in GRPS if len(VAL[lab][g])}
    out = []
    for (g1, g2), p in pairs:
        lo, hi = (g1, g2) if med[g1] <= med[g2] else (g2, g1)
        out.append((med[lo], med[hi], lo, hi, p))
    return sorted(out, key=lambda t: t[1] - t[0])


fig = plt.figure(figsize=(7.3, 3.35), facecolor="white")
gs = fig.add_gridspec(1, 2, wspace=0.20, left=0.072, right=0.995, top=0.735,
                      bottom=0.215, width_ratios=[1.18, 1.0])
fig.text(0.5, 0.982, "Distance from SPP1+ macrophages to exhausted T cells",
         ha="center", va="top", fontsize=8.6)
fig.text(0.5, 0.912, "each point = one section (median over all SPP1+ Mac in that section)",
         ha="center", va="top", fontsize=6.2, color="0.30")
for j, g in enumerate(GRPS):
    fig.add_artist(plt.Line2D([0.285 + j * 0.165], [0.838], marker="s", ls="", ms=4.6,
                              mfc=GCOL[g], mec="black", mew=0.4, transform=fig.transFigure))
    fig.text(0.307 + j * 0.165, 0.838, GLAB[g], fontsize=6.4, va="center", ha="left")

x = np.arange(len(SUBSETS))
W = 0.25
offs = {g: (i - 1) * W for i, g in enumerate(GRPS)}

# ---------------- Panel A: grouped boxplots on a shared axis ----------------
axA = fig.add_subplot(gs[0, 0])
allv = [VAL[lab][g] for _, lab in SUBSETS for g in GRPS if len(VAL[lab][g])]
geo = fs.tight_ylim(axA, allv, n_brackets=3, show_clip_marker=False)
for i, (key, lab) in enumerate(SUBSETS):
    d, kw, pairs = stats(lab)
    for g in GRPS:
        v = VAL[lab][g]
        if len(v) == 0:
            continue
        bp = axA.boxplot([v], positions=[i + offs[g]], widths=W * 0.82, patch_artist=True,
                         showfliers=False,
                         medianprops=dict(color="black", linewidth=0.9),
                         boxprops=dict(linewidth=0.45, edgecolor="black"),
                         whiskerprops=dict(linewidth=0.45, color="black"),
                         capprops=dict(linewidth=0.45, color="black"))
        bp["boxes"][0].set_facecolor(GCOL[g])
        bp["boxes"][0].set_alpha(0.85)
        rng = np.random.default_rng(42)
        axA.scatter(i + offs[g] + rng.normal(0, 0.028, len(v)), v, s=3.6, color="black",
                    alpha=0.5, linewidths=0, zorder=3)
        over = v[v > geo["top"]]                 # points above the axis top: filled triangle at the top
        if len(over):
            axA.scatter(i + offs[g] + np.linspace(-0.035, 0.035, len(over)),
                        np.full(len(over), geo["top"]), marker="^", s=13, color="black",
                        linewidths=0, clip_on=False, zorder=8)
    axA.text(i, -0.115, f"n = {len(d.get('Stable', []))} / {len(d.get('Growth', []))} / "
                        f"{len(d.get('Solid', []))}", transform=axA.get_xaxis_transform(),
             ha="center", va="top", fontsize=5.2, color="0.40")
    by_pair = {tuple(sorted([g1, g2])): p for (g1, g2), p in pairs}
    order = sorted(itertools.combinations(GRPS, 2),
                   key=lambda c: abs(offs[c[1]] - offs[c[0]]))   # narrow below, wide above
    for k, (g1, g2) in enumerate(order):
        fs._draw_bracket(axA, i + offs[g1], i + offs[g2],
                         geo["base"] + k * geo["step_y"],
                         fmt(by_pair[tuple(sorted([g1, g2]))]), None, 0.026, 5.2, 0.5)
axA.set_xticks(x)
axA.set_xticklabels([lab for _, lab in SUBSETS])
axA.set_xlim(-0.52, len(SUBSETS) - 0.48)
axA.set_ylabel("Distance to nearest T cell ($\\mu$m)")
axA.set_title("Three subsets on a shared axis", fontsize=7.6, pad=12)
fs.classic(axA)

# ---------------- Panel B: group-median lines + per-subset vertical significance ----------------
axB = fig.add_subplot(gs[0, 1])
med_all = []
for g in GRPS:
    med = [float(np.median(VAL[lab][g])) if len(VAL[lab][g]) else np.nan
           for _, lab in SUBSETS]
    med_all += [v for v in med if v == v]
    axB.plot(x, med, marker="o", ms=3.6, lw=1.6, color=GCOL[g], mfc="white", mew=0.9,
             zorder=5)
axB.set_ylim(0, max(med_all) * 1.30)
for i, (key, lab) in enumerate(SUBSETS):
    d, kw, pairs = stats(lab)
    for k, (y_lo, y_hi, g_lo, g_hi, p) in enumerate(nest(lab, pairs)):
        fs._draw_vbracket(axB, i + 0.20 + k * 0.095, y_lo, y_hi, fmt(p), fs=5.0,
                          lw=0.5, cap=0.020, txt_dx=0.040)
    axB.text(i, -0.115, "KW " + ("P < 0.001" if kw < 0.001 else f"P = {kw:.3f}"),
             transform=axB.get_xaxis_transform(), ha="center", va="top",
             fontsize=5.4, color="0.25")
axB.set_xticks(x)
axB.set_xticklabels([lab for _, lab in SUBSETS])
axB.set_xlim(-0.30, len(SUBSETS) - 0.18)
axB.set_title("Group medians, one line per group", fontsize=7.6, pad=12)
axB.set_ylabel("")
fs.classic(axB)

fig.text(0.5, 0.052,
         "P values: two-sided Mann–Whitney U with Bonferroni correction for 3 comparisons "
         "(vertical markers / brackets);  Kruskal–Wallis test below each subset.\n"
         "n = number of sections; a value is undefined where a section contains no SPP1+ Mac "
         "(0 or 1 cell) and is reported as missing, not as a drop.",
         ha="center", va="bottom", fontsize=5.2, color="0.35", linespacing=1.6)

fig.savefig(os.path.join(OUT, "05_distance_shared_axis.pdf"))
fig.savefig(os.path.join(OUT, "05_distance_shared_axis.png"))
fig.savefig(os.path.join(OUT_FORMAL, f"{NAME_FORMAL}.pdf"))
fig.savefig(os.path.join(OUT_FORMAL, f"{NAME_FORMAL}.png"))
plt.close(fig)

print("Panel A axis range:", {k: (round(v, 1) if isinstance(v, (int, float)) else v)
                              for k, v in geo.items() if k != "over"})
for key, lab in SUBSETS:
    d, kw, pairs = stats(lab)
    print(f"{lab:10s} medians " + "  ".join(
        f"{g}={np.median(VAL[lab][g]):.0f}(n={len(VAL[lab][g])})" for g in GRPS)
        + f"   KW p={kw:.4f}  " + "  ".join(f"{a}-{b}:{p:.4f}" for (a, b), p in pairs))
print("saved:", os.path.join(OUT, "05_distance_shared_axis.png"))
print("saved:", os.path.join(OUT_FORMAL, f"{NAME_FORMAL}.png"))
