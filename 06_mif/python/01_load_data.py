#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# mIF step 1 - load the single-cell data and compute per-section metrics
#
# Purpose   : Load every HALO single-cell coordinate CSV for both antibody panels
#             (P6 and P7), match each section to its growth group, and compute the
#             per-section cell metrics used downstream: marker percentage, marker
#             density, tumour- and stroma-stratified percentages and the composite
#             marker ratios. The raw per-section tables are saved as pickles and
#             the metric tables as CSVs.
#
# Inputs    : <PROJECT_ROOT>/mif_batch1_cells/{region_TC,region_IM,region_tumor}/*.csv
#             <PROJECT_ROOT>/mif_batch2_cells/{TC,IM,TUMOR,CD20+}/*.csv
#             <PROJECT_ROOT>/sample_id_map.xlsx   (sheets batch1_map, batch2_map)
#             vendor Classified / TUMOR / STROMA area sheets
#             (VENDOR_P6_STATS and VENDOR_P7_STATS in config.py)
#
# Outputs   : <OUTPUT_DATA>/P6_{TC,IM,Tumor}_metrics.csv
#             <OUTPUT_DATA>/P7_{TC,IM,TUMOR,CD20+}_metrics.csv
#             <OUTPUT_DATA>/{p6_data,p7_data,p6_metrics,p7_metrics}.pkl
#
# Usage     : PROJECT_ROOT=/path/to/data python 01_load_data.py
# -----------------------------------------------------------------------------

import os, pickle, warnings, re
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')

from config import *

# Pixel size: exported single-cell coordinates are in PIXELS, not microns.
# Each pixel = 0.3455006 μm (confirmed against company NND data, ratio≈0.99).
PX_TO_UM = 0.3455006
np.random.seed(42)  # reproducibility

def find_csv_files(base_dir):
    """Recursively find all CSV files excluding ._ hidden files."""
    csv_files = []
    for root, dirs, files in os.walk(base_dir):
        for f in files:
            if f.endswith('.csv') and not f.startswith('._'):
                csv_files.append(os.path.join(root, f))
    return csv_files

# Filename ID typos in the P7 exports:
#   P7 files say 'BT007050-9' but the sample is BT0007050 (Solid)
#   P7 files say 'S0419341-5' but the sample is S0419314 (Solid)
ID_ALIASES = {'BT007050': 'BT0007050', 'S0419341': 'S0419314'}

def extract_sample_id(filepath):
    """Extract simplified sample ID from filename."""
    fname = os.path.basename(filepath)
    # Try to extract the ID pattern like BT0004455-9 or S0367740-5 or TS0006760-6.
    # NOTE: TS must be tried BEFORE S — the S pattern would otherwise match
    # the 'S' inside 'TS0006760-6' and silently drop the TS samples.
    patterns = [
        r'(BT\d+-\d+)', r'(TS\d+-\d+)', r'(S\d+-\d+)',
        r'(BT\d+)', r'(TS\d+)', r'(S\d+)',
        r'(BT\d+-\w+)', r'(TS\d+-\w+)', r'(S\d+-\w+)',
    ]
    for pat in patterns:
        m = re.search(pat, fname)
        if m:
            sid = m.group(1)
            # normalize company typos in the base ID before matching
            for bad, good in ID_ALIASES.items():
                if sid.startswith(bad):
                    sid = sid.replace(bad, good, 1)
                    break
            return sid
    return None

def load_mapping(sheet_name):
    """Load sample-to-group mapping from the given sheet, fixing 'Satble' typo.

    The workbook has one sheet per panel: 'batch1_map' (31 samples) and
    'batch2_map' (32 samples incl. S0346026, which has P7 data only —
    vendor-confirmed). Each panel must be matched against its own sheet."""
    df = pd.read_excel(MAPPING_FILE, sheet_name=sheet_name)
    df['Group_clean'] = df['Group'].str.strip().replace({'Satble': 'Stable'})
    mapping = {}
    for _, row in df.iterrows():
        sid = str(row['image_tag']).strip()
        mapping[sid] = row['Group_clean']
    return mapping

def match_sample_id(extracted_sid, mapping):
    """Match extracted sample ID to mapping.
    extracted: 'BT0004455-9', mapping key: 'BT0004455' -> startswith match."""
    if extracted_sid in mapping:
        return extracted_sid
    # Try prefix matching
    for map_id in mapping:
        if extracted_sid.startswith(map_id) and len(map_id) >= 6:
            return map_id
        if map_id.startswith(extracted_sid) and len(extracted_sid) >= 6:
            return map_id
    return None

# ============================================================
# Vendor (company) HALO area lookup — density denominators
# ============================================================
# The company's density = cell count / Classified Area (TUMOR+STROMA tissue area,
# excluding blank). This area comes from HALO's tissue classifier and CANNOT be
# recomputed from the cell CSVs (no tissue polygon exported), so it is read from
# the vendor's official statistics xlsx.

VENDOR_AREA_ROWS = {'classified': 'Classified Area (mm²)',
                    'tumor': 'TUMOR Area (mm²)',
                    'stroma': 'STROMA Area (mm²)'}

# vendor sheet name -> our region key (P6 uses 'Tumor', P7 uses 'TUMOR')
# NOTE: keys are the SHEET names in the vendor xlsx, values are our region keys.
P6_VENDOR_SHEETS = {'TC': 'TC', 'IM': 'IM', 'TUMOR': 'Tumor'}
P7_VENDOR_SHEETS = {'TC': 'TC', 'IM': 'IM', 'TUMOR': 'TUMOR', 'CD20+': 'CD20+'}

def _base_id(sid):
    """Strip the '-N' suffix: 'TS0006760-6' -> 'TS0006760', 'S0444530-7' -> 'S0444530'."""
    return re.sub(r'-\d+$', '', str(sid).strip())

def load_vendor_areas(stats_path, sheet_to_region):
    """Return {region_key: {base_id: {'classified':..,'tumor':..,'stroma':..}}}."""
    xl = pd.ExcelFile(stats_path)
    areas = {}
    for sheet, region_key in sheet_to_region.items():
        if sheet not in xl.sheet_names:
            print(f"  WARNING: vendor sheet [{sheet}] not found in {os.path.basename(stats_path)}")
            continue
        df = pd.read_excel(stats_path, sheet_name=sheet)
        labels = df['Image Tag'].astype(str).str.strip()
        def get_row(label):
            hits = labels[labels == label].index
            return hits[0] if len(hits) else None
        r_class, r_tumor, r_stroma = (get_row(VENDOR_AREA_ROWS[k]) for k in ('classified', 'tumor', 'stroma'))
        reg = {}
        for col in df.columns[1:]:
            m = re.search(r'-((?:BT|TS|S)\d+)(?:-\d+)?-MERGE', col, re.IGNORECASE)
            if not m:
                continue
            base = m.group(1)
            base = ID_ALIASES.get(base, base)
            try:
                reg[base] = {
                    'classified': float(df.iloc[r_class][col]),
                    'tumor': float(df.iloc[r_tumor][col]),
                    'stroma': float(df.iloc[r_stroma][col]),
                }
            except (TypeError, ValueError):
                continue
        areas[region_key] = reg
        print(f"  Vendor areas loaded: sheet[{sheet}] -> region[{region_key}]: {len(reg)} samples")
    return areas

def load_all_p6_data(mapping):
    """Load all P6 CSV files from TC, IM, Tumor directories."""
    all_data = {'TC': {}, 'IM': {}, 'Tumor': {}}

    region_map = {
        'region_TC': 'TC',
        'region_IM': 'IM',
        'region_tumor': 'Tumor'
    }

    for dir_name, region_key in region_map.items():
        dir_path = os.path.join(P6_BASE, dir_name)
        csv_files = find_csv_files(dir_path)
        print(f"P6-{region_key}: found {len(csv_files)} CSV files")
        unmatched = []

        for fp in csv_files:
            sid = extract_sample_id(fp)
            if sid is None:
                unmatched.append((fp, 'no ID pattern'))
                continue
            map_id = match_sample_id(sid, mapping)
            if map_id is None:
                unmatched.append((fp, f'ID {sid} not in mapping'))
                continue
            try:
                df = pd.read_csv(fp, encoding='utf-8')
            except:
                try:
                    df = pd.read_csv(fp, encoding='latin-1')
                except:
                    print(f"  SKIP encoding error: {os.path.basename(fp)[:60]}")
                    continue
            all_data[region_key][map_id] = {
                'group': mapping[map_id],
                'df': df,
                'file': fp
            }

        if unmatched:
            print(f"  P6-{region_key}: {len(unmatched)} UNMATCHED files (dropped):")
            for fp, why in unmatched:
                print(f"    - {os.path.basename(fp)[:70]} ({why})")

    for region_key in all_data:
        n_stable = sum(1 for v in all_data[region_key].values() if v['group']=='Stable')
        n_growth = sum(1 for v in all_data[region_key].values() if v['group']=='Growth')
        n_solid = sum(1 for v in all_data[region_key].values() if v['group']=='Solid')
        print(f"  P6-{region_key}: Stable={n_stable}, Growth={n_growth}, Solid={n_solid}")

    return all_data

def load_all_p7_data(mapping):
    """Load all P7 CSV files."""
    all_data = {'TC': {}, 'IM': {}, 'TUMOR': {}, 'CD20+': {}}

    region_dirs = {
        'TC': os.path.join(P7_BASE, 'TC'),
        'IM': os.path.join(P7_BASE, 'IM'),
        'TUMOR': os.path.join(P7_BASE, 'TUMOR'),
        'CD20+': os.path.join(P7_BASE, 'CD20+')
    }

    for region_key, dir_path in region_dirs.items():
        csv_files = find_csv_files(dir_path)
        print(f"P7-{region_key}: found {len(csv_files)} CSV files")
        unmatched = []

        for fp in csv_files:
            sid = extract_sample_id(fp)
            if sid is None:
                unmatched.append((fp, 'no ID pattern'))
                continue
            map_id = match_sample_id(sid, mapping)
            if map_id is None:
                unmatched.append((fp, f'ID {sid} not in mapping'))
                continue
            try:
                df = pd.read_csv(fp, encoding='utf-8')
            except:
                try:
                    df = pd.read_csv(fp, encoding='latin-1')
                except:
                    continue
            all_data[region_key][map_id] = {
                'group': mapping[map_id],
                'df': df,
                'file': fp
            }

        if unmatched:
            print(f"  P7-{region_key}: {len(unmatched)} UNMATCHED files (dropped):")
            for fp, why in unmatched:
                print(f"    - {os.path.basename(fp)[:70]} ({why})")

    for region_key in all_data:
        n_stable = sum(1 for v in all_data[region_key].values() if v['group']=='Stable')
        n_growth = sum(1 for v in all_data[region_key].values() if v['group']=='Growth')
        n_solid = sum(1 for v in all_data[region_key].values() if v['group']=='Solid')
        print(f"  P7-{region_key}: Stable={n_stable}, Growth={n_growth}, Solid={n_solid}")

    return all_data

def compute_p6_metrics(sample_data, sample_id, region_key, area_map=None):
    """Compute per-sample metrics for P6 panel."""
    df = sample_data['df']
    total = len(df)

    # Density denominator = vendor Classified Area (TUMOR+STROMA, no blank).
    # Fallback: bounding-box area (px->μm->mm²) if vendor area is missing.
    x_range = df['XMax'].max() - df['XMin'].min()
    y_range = df['YMax'].max() - df['YMin'].min()
    bbox_mm2 = (x_range * y_range) * (PX_TO_UM ** 2) / 1e6 if x_range > 0 and y_range > 0 else 1.0
    vend = area_map.get(_base_id(sample_id), {}) if area_map else {}
    if vend and vend.get('classified'):
        area_mm2 = vend['classified']
    else:
        area_mm2 = bbox_mm2
        print(f"  WARNING: no vendor area for {region_key}/{sample_id}; fallback bbox={bbox_mm2:.3f} mm²")

    metrics = {
        'sample_id': sample_id,
        'group': sample_data['group'],
        'region': region_key,
        'total_cells': total,
        'area_mm2': area_mm2,
    }

    # Percentage and density for each marker
    for col in P6_ALL_MARKER_COLS:
        if col in df.columns:
            count = (df[col] == 1).sum()
            pct = count / total * 100
            density = count / area_mm2
            metrics[f"{col}_pct"] = pct
            metrics[f"{col}_density"] = density

    # Composite immune score
    cd8_pct = metrics.get('CD8 Positive Classification_pct', 0)
    cd4_pct = metrics.get('CD4 Positive Classification_pct', 0)
    cd20_pct = metrics.get('CD20 Positive Classification_pct', 0)
    cd68_pct = metrics.get('CD68 Positive Classification_pct', 0)
    metrics['Immune_Score'] = (cd8_pct + cd4_pct + cd20_pct + cd68_pct) / 4

    # GZMB/PD1 ratio
    gzmb_cd8 = metrics.get('GZMB+CD8+_pct', 0)
    pd1_cd8 = metrics.get('PD1+CD8+_pct', 0)
    metrics['GZMB_PD1_CD8_ratio'] = gzmb_cd8 / pd1_cd8 if pd1_cd8 > 0 else np.nan

    # Tumor purity
    metrics['Tumor_Purity'] = metrics.get('PANCK Positive Classification_pct', 0)
    metrics['Stroma_Score'] = metrics.get('aSMA Positive Classification_pct', 0)

    # Subset proportions
    cd8_total = metrics.get('CD8 Positive Classification_pct', 0)
    cd4_total = metrics.get('CD4 Positive Classification_pct', 0)
    metrics['GZMB_CD8_of_CD8'] = (gzmb_cd8 / cd8_total * 100) if cd8_total > 0 else 0
    metrics['PD1_CD8_of_CD8'] = (pd1_cd8 / cd8_total * 100) if cd8_total > 0 else 0
    metrics['CD8_CD4_ratio'] = (cd8_total / cd4_total) if cd4_total > 0 else np.nan

    gzmb_cd4 = metrics.get('GZMB+CD4+_pct', 0)
    pd1_cd4 = metrics.get('PD1+CD4+_pct', 0)
    metrics['GZMB_CD4_of_CD4'] = (gzmb_cd4 / cd4_total * 100) if cd4_total > 0 else 0
    metrics['PD1_CD4_of_CD4'] = (pd1_cd4 / cd4_total * 100) if cd4_total > 0 else 0

    # SPP1/CD68 ratio
    spp1_cd68 = metrics.get('SPP1+CD68+_pct', 0)
    cd68_total_pct = metrics.get('CD68 Positive Classification_pct', 0)
    metrics['SPP1_CD68_ratio'] = (spp1_cd68 / cd68_total_pct * 100) if cd68_total_pct > 0 else 0

    # B cell fraction of immune
    immune_sum = cd8_total + cd4_total + cd20_pct + cd68_total_pct
    metrics['B_cell_immune_fraction'] = (cd20_pct / immune_sum * 100) if immune_sum > 0 else 0

    # TUMOR/STROMA stratified (percentages within compartment + compartment density)
    if LABEL_COL in df.columns:
        for compartment in ['TUMOR', 'STROMA']:
            comp_df = df[df[LABEL_COL] == compartment]
            comp_total = len(comp_df)
            comp_area = vend.get(compartment.lower()) if vend else None
            if comp_total > 0:
                for col in P6_ALL_MARKER_COLS:
                    if col in df.columns:
                        count = (comp_df[col] == 1).sum()
                        pct = count / comp_total * 100
                        metrics[f"{col}_{compartment}_pct"] = pct
                        if comp_area:
                            metrics[f"{col}_{compartment}_density"] = count / comp_area

    return metrics

def compute_p7_metrics(sample_data, sample_id, region_key, area_map=None):
    """Compute per-sample metrics for P7 panel."""
    df = sample_data['df']
    total = len(df)

    x_range = df['XMax'].max() - df['XMin'].min()
    y_range = df['YMax'].max() - df['YMin'].min()
    bbox_mm2 = (x_range * y_range) * (PX_TO_UM ** 2) / 1e6 if x_range > 0 and y_range > 0 else 1.0
    vend = area_map.get(_base_id(sample_id), {}) if area_map else {}
    if vend and vend.get('classified'):
        area_mm2 = vend['classified']
    else:
        area_mm2 = bbox_mm2
        print(f"  WARNING: no vendor area for {region_key}/{sample_id}; fallback bbox={bbox_mm2:.3f} mm²")

    metrics = {
        'sample_id': sample_id,
        'group': sample_data['group'],
        'region': region_key,
        'total_cells': total,
        'area_mm2': area_mm2,
    }

    for col in P7_ALL_MARKER_COLS:
        if col in df.columns:
            count = (df[col] == 1).sum()
            pct = count / total * 100
            density = count / area_mm2
            metrics[f"{col}_pct"] = pct
            metrics[f"{col}_density"] = density

    # HLA-DR/CD68 ratio
    hladr_cd68 = metrics.get('HLA-DR+CD68+_pct', 0)
    cd68_pct = metrics.get('CD68 Positive Classification_pct', 0)
    metrics['HLADR_CD68_ratio'] = (hladr_cd68 / cd68_pct * 100) if cd68_pct > 0 else 0

    # SPP1/CCL19 ratio (need cross-panel later)

    if LABEL_COL in df.columns:
        for compartment in ['TUMOR', 'STROMA']:
            comp_df = df[df[LABEL_COL] == compartment]
            comp_total = len(comp_df)
            comp_area = vend.get(compartment.lower()) if vend else None
            if comp_total > 0:
                for col in P7_ALL_MARKER_COLS:
                    if col in df.columns:
                        count = (comp_df[col] == 1).sum()
                        pct = count / comp_total * 100
                        metrics[f"{col}_{compartment}_pct"] = pct
                        if comp_area:
                            metrics[f"{col}_{compartment}_density"] = count / comp_area

    return metrics

# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    print("=" * 60)
    print("Phase 1: Loading all data and computing metrics")
    print("=" * 60)

    # Load panel-specific mappings. S0346026 (Solid) appears only in the P7
    # sheet: it has P7 data but no P6 data (vendor-confirmed), so the P6
    # loader must NOT see it — otherwise the P6 panel would silently drop it
    # as unmatched, and the P7 panel would miss it entirely.
    mapping_p6 = load_mapping('batch1_map')
    print(f"Mapping (P6 sheet) loaded: {len(mapping_p6)} samples")
    print(f"Groups: {pd.Series(list(mapping_p6.values())).value_counts().to_dict()}")

    # Load P6
    print("\n--- P6 ---")
    p6_data = load_all_p6_data(mapping_p6)

    mapping_p7 = load_mapping('batch2_map')
    print(f"Mapping (P7 sheet) loaded: {len(mapping_p7)} samples")
    print(f"Groups: {pd.Series(list(mapping_p7.values())).value_counts().to_dict()}")

    # Load P7
    print("\n--- P7 ---")
    p7_data = load_all_p7_data(mapping_p7)

    # Load vendor (company) Classified/TUMOR/STROMA areas — density denominators
    print("\n--- Loading vendor areas ---")
    p6_areas = load_vendor_areas(VENDOR_P6_STATS, P6_VENDOR_SHEETS)
    p7_areas = load_vendor_areas(VENDOR_P7_STATS, P7_VENDOR_SHEETS)

    # Compute metrics
    print("\n--- Computing P6 metrics ---")
    p6_metrics = {'TC': [], 'IM': [], 'Tumor': []}
    for region in ['TC', 'IM', 'Tumor']:
        for sid, sdata in p6_data[region].items():
            m = compute_p6_metrics(sdata, sid, region, area_map=p6_areas.get(region, {}))
            p6_metrics[region].append(m)
        df = pd.DataFrame(p6_metrics[region])
        print(f"  P6-{region}: {len(df)} samples, {len(df.columns)} metrics")

    print("\n--- Computing P7 metrics ---")
    p7_metrics = {'TC': [], 'IM': [], 'TUMOR': [], 'CD20+': []}
    for region in ['TC', 'IM', 'TUMOR', 'CD20+']:
        for sid, sdata in p7_data[region].items():
            m = compute_p7_metrics(sdata, sid, region, area_map=p7_areas.get(region, {}))
            p7_metrics[region].append(m)
        df = pd.DataFrame(p7_metrics[region])
        print(f"  P7-{region}: {len(df)} samples, {len(df.columns)} metrics")

    # Save intermediate data
    os.makedirs(OUTPUT_DATA, exist_ok=True)

    for region in ['TC', 'IM', 'Tumor']:
        df = pd.DataFrame(p6_metrics[region])
        df.to_csv(os.path.join(OUTPUT_DATA, f"P6_{region}_metrics.csv"), index=False)

    for region in ['TC', 'IM', 'TUMOR', 'CD20+']:
        df = pd.DataFrame(p7_metrics[region])
        df.to_csv(os.path.join(OUTPUT_DATA, f"P7_{region}_metrics.csv"), index=False)

    # Save raw data references
    with open(os.path.join(OUTPUT_DATA, 'p6_data.pkl'), 'wb') as f:
        pickle.dump(p6_data, f)
    with open(os.path.join(OUTPUT_DATA, 'p7_data.pkl'), 'wb') as f:
        pickle.dump(p7_data, f)
    with open(os.path.join(OUTPUT_DATA, 'p6_metrics.pkl'), 'wb') as f:
        pickle.dump(p6_metrics, f)
    with open(os.path.join(OUTPUT_DATA, 'p7_metrics.pkl'), 'wb') as f:
        pickle.dump(p7_metrics, f)

    print("\nPhase 1 complete. Data saved to", OUTPUT_DATA)
