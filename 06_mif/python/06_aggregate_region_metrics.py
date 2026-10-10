#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF step 6 - per-section metrics of the CD20+ aggregate region
#
# Purpose   : Divide each per-section measurement of the annotated CD20+ B-cell
#             aggregate region by the section's circled aggregate count, and report
#             the group comparisons.
#
# Inputs    : <VENDOR_DIR>/vendor_cd20_aggregate_counts_batch2.xlsx
#             <VENDOR_DIR>/vendor_cell_stats_batch2.xlsx, sheet 'CD20+'
#             rows used: 'Classified Area (mm2)', 'Total Cells',
#             'CD20 Positive Cells', 'CCL19 Positive Cells'
#             <PROJECT_ROOT>/sample_id_map.xlsx, sheet 'batch2_map'
#
# Outputs   : <OUTPUT_DATA>/aggregate_region_metrics_CD20.csv
#             one row per section: sample_id, group, n_aggregates, ccl19_per_agg
#             and the vendor per-section measurements the ratio is built from
#
# Usage     : PROJECT_ROOT=/path/to/data python 06_aggregate_region_metrics.py
# -----------------------------------------------------------------------------

import os, re, warnings
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
warnings.filterwarnings('ignore')

GROUP_ORDER = ['Stable', 'Growth', 'Solid']
GROUP_LABELS = {'Stable': 'Stable', 'Growth': 'GGO-growth', 'Solid': 'SOLID-growth'}

PROJECT_ROOT = os.environ.get("PROJECT_ROOT", "/path/to/GGN_progression_data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "mIF_Validation_Results")
DATA_DIR = os.path.join(OUTPUT_DIR, "data")
VENDOR_DIR = os.path.join(OUTPUT_DIR, "vendor_reference")
MAPPING_FILE = os.path.join(PROJECT_ROOT, "sample_id_map.xlsx")
VENDOR_AGG = os.path.join(VENDOR_DIR, "vendor_cd20_aggregate_counts_batch2.xlsx")
VENDOR_STATS = os.path.join(VENDOR_DIR, "vendor_cell_stats_batch2.xlsx")

# Two sections are spelled differently in the two vendor exports; both
# spellings denote the same block.
ID_ALIASES = {'BT007050': 'BT0007050', 'S0419341': 'S0419314'}

# Rows of the 'CD20+' sheet, by their label in the first column.
ROW_AREA = 'Classified Area'
ROW_TOTAL = 'Total Cells'
ROW_CD20 = 'CD20 Positive Cells'
ROW_CCL19 = 'CCL19 Positive Cells'


def tag_section(tag):
    """('S0462313', 6) from an image tag such as ...-S0462313-6-MERGE.ome.tiff."""
    m = re.search(r'-((?:BT|TS|S)\d+)-(\d+)-', str(tag), re.I)
    if not m:
        return None, None
    base = m.group(1).upper()
    return ID_ALIASES.get(base, base), int(m.group(2))


def section_id(tag):
    base, sec = tag_section(tag)
    return None if base is None else f"{base}-{sec}"


def find_row(df, label):
    """First row whose label starts with `label`, compared after case, spacing
    and punctuation are removed ('CD20 Positive Cells' -> 'cd20positivecells').

    Digits are part of the key: the sheet lists CD3, CD20, CD31, CD68 and CA9
    side by side, so dropping them would make 'CD20 Positive Cells' match
    'CD68 Positive Cells' instead. The 'Analysis Region' block is the first of
    the sheet's several blocks, so the first match is the whole-region figure;
    rows beginning with '%' are the percentage twins and are skipped.
    """
    want = re.sub(r'[^a-z0-9]', '', label.lower())
    for i in range(len(df)):
        raw = str(df.iloc[i, 0])
        if raw.lstrip().startswith('%'):
            continue
        if re.sub(r'[^a-z0-9]', '', raw.lower()).startswith(want):
            return i
    raise KeyError(f"row '{label}' not found in the CD20+ sheet")


def get_group_data(df, col):
    return {g: df.loc[df['group'] == g, col].dropna().values.astype(float)
            for g in GROUP_ORDER}


def report(df, col):
    """Kruskal-Wallis plus the three pairwise two-sided Mann-Whitney U tests,
    Bonferroni-corrected across the comparisons actually made."""
    dd = get_group_data(df, col)
    present = [g for g in GROUP_ORDER if len(dd[g]) >= 2]
    if len(present) < 2:
        return
    print(f"  {col}: KW P = {kruskal(*[dd[g] for g in present]).pvalue:.4g}")
    pairs = [(present[i], present[j])
             for i in range(len(present)) for j in range(i + 1, len(present))]
    for g1, g2 in pairs:
        p = mannwhitneyu(dd[g1], dd[g2], alternative='two-sided').pvalue
        print(f"    {GROUP_LABELS[g1]} vs {GROUP_LABELS[g2]}: "
              f"P = {min(p * len(pairs), 1.0):.4g} (Bonferroni x{len(pairs)}; "
              f"uncorrected {p:.4g}, n = {len(dd[g1])} vs {len(dd[g2])})")


# ============================================================
# 1. Circled-aggregate count per section (vendor)
# ============================================================
print("Reading vendor circled-aggregate counts ...", flush=True)
vagg = pd.read_excel(VENDOR_AGG)
vagg.columns = ['image_tag', 'n_aggregates']
vagg['sample_id'] = vagg['image_tag'].apply(section_id)
vagg = vagg.dropna(subset=['sample_id'])
counts = vagg.groupby('sample_id')['n_aggregates'].first()
print(f"  {len(counts)} sections with a circled-aggregate count", flush=True)

# ============================================================
# 2. Per-section measurements of the CD20+ region (vendor)
# ============================================================
print("Reading per-section measurements of the CD20+ region ...", flush=True)
sdf = pd.read_excel(VENDOR_STATS, sheet_name='CD20+')
i_area, i_total = find_row(sdf, ROW_AREA), find_row(sdf, ROW_TOTAL)
i_cd20, i_ccl19 = find_row(sdf, ROW_CD20), find_row(sdf, ROW_CCL19)
print(f"  rows used: area={i_area}, total={i_total}, CD20={i_cd20}, CCL19={i_ccl19}", flush=True)

rows = []
for col in sdf.columns[1:]:
    sid = section_id(col)
    if sid is None:
        continue
    rows.append({
        'sample_id': sid,
        'region_area_mm2': float(sdf.iloc[i_area][col]),
        'region_total_cells': float(sdf.iloc[i_total][col]),
        'region_cd20_cells': float(sdf.iloc[i_cd20][col]),
        'region_ccl19_cells': float(sdf.iloc[i_ccl19][col]),
    })
region = pd.DataFrame(rows).drop_duplicates('sample_id', keep='first')
print(f"  {len(region)} sections in the CD20+ region sheet", flush=True)

# ============================================================
# 3. Groups, then the per-aggregate division
# ============================================================
map_df = pd.read_excel(MAPPING_FILE, sheet_name='batch2_map')
map_df['Group_clean'] = map_df['Group'].str.strip().replace({'Satble': 'Stable'})
# The mapping workbook keys a section either as the bare block code or with its
# section suffix; the group belongs to the block, so both sides are reduced to
# the base code before they are joined.
map_df['base_id'] = (map_df['image_tag'].astype(str).str.strip()
                     .str.replace(r'-\d+$', '', regex=True))
id2g = dict(zip(map_df['base_id'], map_df['Group_clean']))

metrics = region.merge(counts.rename('n_aggregates'), on='sample_id', how='left')
metrics['base_id'] = metrics['sample_id'].str.replace(r'-\d+$', '', regex=True)
metrics['group'] = metrics['base_id'].map(id2g)
unmapped = metrics[metrics['group'].isna()]
if len(unmapped):
    print("  Sections dropped (absent from the mapping, never silently): "
          + ", ".join(unmapped['sample_id']), flush=True)
metrics = metrics.dropna(subset=['group']).drop(columns=['base_id']).copy()

# A section with no circled aggregate has n_aggregates = 0 (nothing was
# circled); its per-aggregate ratio is undefined and stays empty.
metrics['n_aggregates'] = metrics['n_aggregates'].fillna(0).astype(int)
metrics['ccl19_per_agg'] = (metrics['region_ccl19_cells']
                            / metrics['n_aggregates'].replace(0, np.nan))

metrics = metrics[['sample_id', 'group', 'n_aggregates', 'ccl19_per_agg',
                   'region_area_mm2', 'region_total_cells',
                   'region_cd20_cells', 'region_ccl19_cells']]
metrics.to_csv(os.path.join(DATA_DIR, 'aggregate_region_metrics_CD20.csv'), index=False)

print(f"\nSections analysed: {len(metrics)}", flush=True)
print("Per group: " + ", ".join(
    f"{g}={int((metrics['group'] == g).sum())}" for g in GROUP_ORDER), flush=True)
print("Sections with no circled aggregate: "
      f"{int((metrics['n_aggregates'] == 0).sum())}", flush=True)

# ============================================================
# 4. Group comparisons
# ============================================================
print("\n" + "=" * 66)
print("Median per group:")
for col in ['n_aggregates', 'ccl19_per_agg']:
    meds = {g: metrics[metrics['group'] == g][col].median() for g in GROUP_ORDER}
    print(f"  {col:15s} " + "  ".join(f"{GROUP_LABELS[g]}: {meds[g]:.3g}"
                                       for g in GROUP_ORDER))
print("\nGroup comparisons (two-sided Mann-Whitney U, Bonferroni-corrected):")
for col in ['n_aggregates', 'ccl19_per_agg']:
    report(metrics, col)
print("=" * 66)
