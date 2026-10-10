#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF step 9 - CCL19 / TLS axis
#
# Purpose   : Three-panel CCL19 / TLS axis: CCL19+ cells per CD20+ aggregate, the
#             distance from CCR7+ to the nearest CCL19+ cell, and the CD20+
#             aggregate count per section.
#
# Inputs    : <OUTPUT_DATA>/aggregate_region_metrics_CD20.csv           (panel A)
#             <OUTPUT_DATA>/spatial_results/CCR7_CCL19_NND_median.csv  (panel B)
#             <VENDOR_DIR>/vendor_cd20_aggregate_counts_batch2.xlsx     (panel C)
#             <PROJECT_ROOT>/sample_id_map.xlsx, sheet 'batch2_map'
#             <OUTPUT_DATA>/P7_TC_metrics.csv   (cohort section list)
#
# Outputs   : <FIGS_DIR>/preview/07_ccl19_tls_axis.{pdf,png}
#             <FIGS_DIR>/09_CD20_Aggregates/Fig9_CCL19_Aggregate_Axis.{pdf,png}
#             <OUTPUT_DATA>/Fig9_CCL19_aggregate_axis_stats.csv   (all tests)
#
# Usage     : PROJECT_ROOT=/path/to/data python 09_ccl19_tls_axis.py
# -----------------------------------------------------------------------------

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import figstyle as fs

PROJECT_ROOT = os.environ.get("PROJECT_ROOT", "/path/to/GGN_progression_data")
BASE = os.path.join(PROJECT_ROOT, "mIF_Validation_Results")
DATA = os.path.join(BASE, "data")
VENDOR_DIR = os.path.join(BASE, "vendor_reference")
MAPPING_FILE = os.path.join(PROJECT_ROOT, "sample_id_map.xlsx")
OUT_PREVIEW = os.path.join(BASE, "figures", "preview")
OUT_FORMAL = os.path.join(BASE, "figures", "09_CD20_Aggregates")
NAME_FORMAL = "Fig9_CCL19_Aggregate_Axis"
fs.prepare()

GRPS = ["Stable", "Growth", "Solid"]
GLAB, GCOL = fs.GROUP_LABELS, fs.GROUP_COLORS
ID_ALIASES = {'BT007050': 'BT0007050', 'S0419341': 'S0419314'}


# ---------------- Panel A / B data ----------------
def col_by_group(df, col):
    return {g: df.loc[df["group"] == g, col].dropna().values.astype(float) for g in GRPS}


agg_sum = pd.read_csv(os.path.join(DATA, "aggregate_region_metrics_CD20.csv"))
A = col_by_group(agg_sum, "ccl19_per_agg")

nnd = pd.read_csv(os.path.join(DATA, "spatial_results", "CCR7_CCL19_NND_median.csv"))
B = col_by_group(nnd, "obs_nnd")


# ---------------- Panel C: vendor aggregate count ----------------
def extract_base_id(tag):
    m = re.search(r"-((?:BT|TS|S)\d+)(?:-\d+)?-M", str(tag), re.I)
    if not m:
        return None
    b = m.group(1)
    return ID_ALIASES.get(b, b)


agg_raw = pd.read_excel(os.path.join(VENDOR_DIR, "vendor_cd20_aggregate_counts_batch2.xlsx"))
agg_raw.columns = ["ImageTag", "Count"]

map_df = pd.read_excel(MAPPING_FILE, sheet_name="batch2_map")
map_df["Group_clean"] = map_df["Group"].str.strip().replace({"Satble": "Stable"})
map_df["id_base"] = (map_df["image_tag"].astype(str).str.strip()
                     .str.replace(r"-\d+$", "", regex=True))
id2g = dict(zip(map_df["id_base"], map_df["Group_clean"]))

agg_raw["base_id"] = agg_raw["ImageTag"].apply(extract_base_id)
agg_raw["group"] = agg_raw["base_id"].map(id2g)
agg_all = agg_raw.dropna(subset=["group"]).copy()

# Sections without aggregates scored 0.
cohort_ids = (pd.read_csv(os.path.join(DATA, "P7_TC_metrics.csv"))["sample_id"]
              .astype(str).str.replace(r"-\d+$", "", regex=True))
_have = set(agg_all["base_id"])
zero_added = []
for sid in cohort_ids:
    if sid in _have:
        continue
    g = id2g.get(sid)
    if g is None:
        print(f"    WARNING: {sid} has P7 data but is not in the mapping — group not invented, skipped")
        continue
    agg_all = pd.concat([agg_all, pd.DataFrame([{"base_id": sid, "group": g, "Count": 0}])],
                        ignore_index=True)
    zero_added.append((sid, g))
print(f"Sections scored 0: {len(zero_added)} " + ", ".join(f"{s}({g})" for s, g in zero_added))
C = col_by_group(agg_all, "Count")


# ---------------- Plotting ----------------
PANELS = [
    (A, "CCL19$^+$ cells per\nCD20$^+$ aggregate", "Cells per aggregate"),
    (B, "Distance from CCR7$^+$ to\nnearest CCL19$^+$ cell", "Distance ($\\mu$m)"),
    (C, "CD20$^+$ aggregates\nper section", "Count"),
]

fig = plt.figure(figsize=(7.3, 3.30), facecolor="white")
gs = fig.add_gridspec(1, 3, wspace=0.34, left=0.062, right=0.995, top=0.705,
                      bottom=0.255)

fig.text(0.5, 0.982, "CCL19 organizer signal in CD20$^+$ aggregates",
         ha="center", va="top", fontsize=8.6)
fig.text(0.5, 0.916, "P7 CD20$^+$ region (A) and tumour core (B);  each point = one section",
         ha="center", va="top", fontsize=6.2, color="0.30")

# Legend: group colour + group names (the three panels have different n, so each
# panel carries its own n below the x-axis).
for j, g in enumerate(GRPS):
    fig.add_artist(plt.Line2D([0.300 + j * 0.170], [0.838], marker="s", ls="", ms=4.6,
                              mfc=GCOL[g], mec="black", mew=0.4, transform=fig.transFigure))
    fig.text(0.322 + j * 0.170, 0.838, GLAB[g], fontsize=6.4, va="center", ha="left")

rows = []
for k, (groups, title, ylab) in enumerate(PANELS):
    ax = fig.add_subplot(gs[0, k])
    present = [g for g in GRPS if len(groups[g]) >= 2]
    st = fs.panel(ax, groups, present, style="box", ylabel=ylab, title=title)
    ax.title.set_multialignment("center")
    # Each panel's own n (the three panels have different section counts).
    ax.text(0.5, -0.175, "n = " + " / ".join(str(len(groups[g])) for g in GRPS),
            transform=ax.transAxes, ha="center", va="top", fontsize=5.2, color="0.40")
    for g in present:
        rows.append(dict(panel="ABC"[k], metric=title.replace("\n", " "), group=g,
                         n=len(groups[g]), median=float(np.median(groups[g])),
                         mean=float(np.mean(groups[g])), kruskal_p=st["kw_p"]))
    for (g1, g2), p in st["pairs"]:
        rows.append(dict(panel="ABC"[k], metric=title.replace("\n", " "), group=f"{g1}_vs_{g2}",
                         n=np.nan, median=np.nan, mean=np.nan, kruskal_p=np.nan,
                         pairwise_p_bonferroni=p))

fig.text(0.5, 0.052,
         "P values: two-sided Mann–Whitney U with Bonferroni correction for 3 comparisons "
         "(brackets);  Kruskal–Wallis in each panel title.\n"
         "A: CCL19$^+$ cells in the annotated aggregate region, divided by that section's "
         "circled aggregate count.  B: per-cell nearest-neighbour distance, median per section.\n"
         "C: vendor official aggregate count; sections with no aggregate are scored 0 "
         "(in-house verified), so all 32 cohort sections are shown.  Panel A is zero-inflated "
         "by biology, not by dropout.",
         ha="center", va="bottom", fontsize=5.2, color="0.35", linespacing=1.6)

pd.DataFrame(rows).to_csv(os.path.join(DATA, "Fig9_CCL19_aggregate_axis_stats.csv"), index=False)
for dest, name in [(OUT_PREVIEW, "07_ccl19_tls_axis"), (OUT_FORMAL, NAME_FORMAL)]:
    fig.savefig(os.path.join(dest, f"{name}.pdf"))
    fig.savefig(os.path.join(dest, f"{name}.png"))
plt.close(fig)

for k, (groups, title, _y) in enumerate(PANELS):
    present = [g for g in GRPS if len(groups[g]) >= 2]
    st = fs.group_stats(groups, present)
    print(f"Panel {'ABC'[k]}  " + "  ".join(
        f"{g}: {np.median(groups[g]):.2f} (n={len(groups[g])})" for g in GRPS)
        + f"   KW p={st['kw_p']:.4f}  "
        + "  ".join(f"{a}-{b}:{p:.4f}" for (a, b), p in st["pairs"]))
print("saved:", os.path.join(OUT_PREVIEW, "07_ccl19_tls_axis.png"))
print("saved:", os.path.join(OUT_FORMAL, f"{NAME_FORMAL}.png"))
print("saved:", os.path.join(DATA, "Fig9_CCL19_aggregate_axis_stats.csv"))
