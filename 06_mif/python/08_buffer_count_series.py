#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF step 8 - radius series of the PD-1+ T-cell count around SPP1+ macrophages
#
# Purpose   : One panel per T-cell subset (PD-1+CD4+, PD-1+CD8+, PD-1+): group
#             mean +/- SEM across the radii 20, 50 and 100, with a vertical
#             significance marker at each radius.
#
# Inputs    : <OUTPUT_DATA>/Spp1Mac_buffer_contact_TC.csv, written by
#             05_spp1_buffer_contact.py; {CD4,CD8,T}_r{20,50,100} is the number of
#             target cells within radius r of any SPP1+ macrophage
#
# Outputs   : <FIGS_DIR>/preview/06_buffer_count_series.{pdf,png}
#             <FIGS_DIR>/10_SPP1Mac_Buffer_Contact/
#                 Fig10_SPP1Mac_ExhaustedT_Radius_series.{pdf,png}
#             <OUTPUT_DATA>/Fig10_buffer_count_series_stats.csv   (all 27 tests)
#
# Usage     : PROJECT_ROOT=/path/to/data python 08_buffer_count_series.py
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
# Formal version of the same figure.
OUT_FORMAL = os.path.join(BASE, "figures", "10_SPP1Mac_Buffer_Contact")
NAME_FORMAL = "Fig10_SPP1Mac_ExhaustedT_Radius_series"
fs.prepare()

GRPS = ["Stable", "Growth", "Solid"]
GLAB, GCOL = fs.GROUP_LABELS, fs.GROUP_COLORS
TARGETS = [("CD4", "PD-1+CD4+ T"), ("CD8", "PD-1+CD8+ T"), ("T", "PD-1+ T")]
RADII = [20, 50, 100]


raw = pd.read_csv(os.path.join(DATA, "Spp1Mac_buffer_contact_TC.csv"))


def series(tname, r):
    return {g: raw.loc[raw["group"] == g, f"{tname}_r{r}"].dropna().values.astype(float)
            for g in GRPS}


def stats(tname, r):
    from scipy.stats import kruskal, mannwhitneyu
    d = series(tname, r)
    try:
        kw = kruskal(*[v for v in d.values() if len(v) >= 2])[1]
    except ValueError:                     # all zeros (no cell hit): zero variance
        kw = 1.0
    pairs = []
    for g1, g2 in itertools.combinations(GRPS, 2):
        try:
            p = mannwhitneyu(d[g1], d[g2], alternative="two-sided")[1]
        except ValueError:                 # identical values in both groups: scipy errors -> record p = 1
            p = 1.0
        pairs.append(((g1, g2), min(p * 3, 1.0)))
    return d, kw, pairs


def fmt(p):
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}"


fig = plt.figure(figsize=(7.3, 3.35), facecolor="white")
gs = fig.add_gridspec(1, 3, wspace=0.26, left=0.062, right=0.995, top=0.725,
                      bottom=0.315)
fig.text(0.5, 0.982, "Exhausted T cells within a given radius of SPP1+ macrophages",
         ha="center", va="top", fontsize=8.6)
fig.text(0.5, 0.918, "mean cell count per section ± SEM, tumour core (TC)",
         ha="center", va="top", fontsize=6.2, color="0.30")
for j, g in enumerate(GRPS):
    fig.add_artist(plt.Line2D([0.265 + j * 0.205], [0.842], marker="o", ls="-", lw=1.5,
                              ms=4.0, color=GCOL[g], mfc="white", mew=0.9,
                              transform=fig.transFigure))
    fig.text(0.289 + j * 0.205, 0.842, f"{GLAB[g]}  (n = {int((raw['group'] == g).sum())})",
             fontsize=6.4, va="center", ha="left")

x = np.arange(len(RADII))
for k, (tname, tlab) in enumerate(TARGETS):
    ax = fig.add_subplot(gs[0, k])
    ymax = 0.0
    for g in GRPS:
        m = [float(np.mean(series(tname, r)[g])) for r in RADII]
        e = [float(np.std(series(tname, r)[g], ddof=1) / np.sqrt(len(series(tname, r)[g])))
             for r in RADII]
        ymax = max(ymax, max(np.add(m, e)))
        ax.plot(x, m, marker="o", ms=4.2, lw=1.6, color=GCOL[g], mfc="white", mew=0.95,
                zorder=5)
        ax.fill_between(x, np.maximum(np.subtract(m, e), 0), np.add(m, e),
                        color=GCOL[g], alpha=0.16, linewidth=0, zorder=2)
    # Lower the axis floor below 0 so the near-zero lines do not sit on the
    # x-axis, while 0 is still kept as a tick.
    ax.set_ylim(-0.14 * ymax, ymax * 1.30)
    ytop = ax.get_ylim()[1]
    ax.axhline(0, color="0.72", ls=(0, (1.6, 1.8)), lw=0.5, zorder=1)
    for i, r in enumerate(RADII):
        d, kw, pairs = stats(tname, r)
        med = {g: float(np.mean(d[g])) for g in GRPS}
        order = sorted(pairs, key=lambda t: abs(med[t[0][0]] - med[t[0][1]]))
        sig = [t for t in order if t[1] < 0.05]     # only significant pairwise comparisons are drawn
        for j, ((g1, g2), p) in enumerate(sig):
            lo_g, hi_g = (g1, g2) if med[g1] <= med[g2] else (g2, g1)
            lo, hi = med[lo_g], med[hi_g]
            if hi - lo < 0.012 * ytop:              # fully overlapping: keep a short visible marker
                hi = lo + 0.012 * ytop
            fs._draw_vbracket(ax, i + 0.28 + j * 0.13, lo, hi, fmt(p), fs=5.1,
                              lw=0.45, cap=0.020, txt_dx=0.030, color="0.22")
        ax.text(i, -0.145, "KW " + ("P < 0.001" if kw < 0.001 else f"P = {kw:.3f}"),
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=5.3, color="0.25")
    ax.set_xticks(x)
    ax.set_xticklabels([str(r) for r in RADII])
    ax.set_xlabel("Radius around SPP1+ Mac ($\\mu$m)", labelpad=13)
    ax.set_xlim(-0.34, len(RADII) - 0.06)
    ax.set_title(tlab, fontsize=7.6, pad=5)
    if k == 0:
        ax.set_ylabel("Cells per section")
    fs.classic(ax)

fig.text(0.5, 0.048,
         "Line: group mean;  band: SEM;  vertical markers: two-sided Mann–Whitney U, "
         "Bonferroni-corrected for 3 comparisons, drawn only when significant.\n"
         "Kruskal–Wallis below each radius.  Counts are zero-inflated: at 20 and 50 "
         "$\\mu$m most SOLID-growth sections and essentially all Stable / GGO-growth "
         "sections contain no such cell,\n"
         "so many comparisons are tied and return P = 1.  All 27 tests, significant or "
         "not, are listed in the accompanying stats CSV.",
         ha="center", va="bottom", fontsize=5.2, color="0.35", linespacing=1.6)

fig.savefig(os.path.join(OUT, "06_buffer_count_series.pdf"))
fig.savefig(os.path.join(OUT, "06_buffer_count_series.png"))
fig.savefig(os.path.join(OUT_FORMAL, f"{NAME_FORMAL}.pdf"))
fig.savefig(os.path.join(OUT_FORMAL, f"{NAME_FORMAL}.png"))
plt.close(fig)

rows = []
for tname, tlab in TARGETS:
    for r in RADII:
        d, kw, pairs = stats(tname, r)
        for (g1, g2), p in pairs:
            rows.append(dict(target=tlab, radius_um=r, group_1=g1, group_2=g2,
                             mean_1=float(np.mean(d[g1])), mean_2=float(np.mean(d[g2])),
                             n_1=len(d[g1]), n_2=len(d[g2]),
                             p_mannwhitney_bonferroni=p, kruskal_p=kw,
                             shown_in_figure=bool(p < 0.05)))
_stats = pd.DataFrame(rows)
_stats.to_csv(os.path.join(OUT, "06_buffer_count_series_stats.csv"), index=False)
_stats.to_csv(os.path.join(DATA, "Fig10_buffer_count_series_stats.csv"), index=False)

for tname, tlab in TARGETS:
    for r in RADII:
        d, kw, pairs = stats(tname, r)
        print(f"{tlab:12s} r={r:3d}  " + "  ".join(
            f"{g}:{np.mean(d[g]):7.2f}(n={len(d[g])})" for g in GRPS)
            + f"  KW={kw:.4f}  " + "  ".join(f"{a}-{b}:{p:.4f}" for (a, b), p in pairs))
print("saved:", os.path.join(OUT, "06_buffer_count_series.png"))
print("saved:", os.path.join(OUT_FORMAL, f"{NAME_FORMAL}.png"))
print("saved:", os.path.join(DATA, "Fig10_buffer_count_series_stats.csv"))
