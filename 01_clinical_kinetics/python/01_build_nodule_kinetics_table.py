# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Build the analysis-ready nodule kinetics table (999 nodules x 45 columns)
#
# Purpose   : Rebuild the kinetics table every downstream panel reads. The four
#             kinetic columns (VDT, MDT total, MDT solid, annualized CTR) and the
#             grouping, pathology, nodule-type and follow-up columns come from
#             Table S1; the remaining quantitative-CT columns (volume, solid
#             volume, CTR components, change pattern) come from the pipeline
#             kinetics table, matched nodule by nodule.
#
# Inputs    : $PROJECT_ROOT/manuscript/Table/Table S1.xlsx
#             $PROJECT_ROOT/omics/supplementary_clinical_analysis/figure_assembly/rerun_999/input/TableS1_curated_binary_nodule_type_999.xlsx
#             $PROJECT_ROOT/omics/supplementary_clinical_analysis/figure_assembly/rerun_999/code/kinetics.csv
#
# Outputs   : ./kinetics.csv (999 rows x 45 columns) and ./build_reconciliation_log.txt
#
# Usage     : LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 PROJECT_ROOT=/path/to/data \
#                 python3 01_build_nodule_kinetics_table.py
# -----------------------------------------------------------------------------

import os
import csv
import collections

import pandas as pd

# ======================================================================
# Path configuration — every input hangs off PROJECT_ROOT
# ======================================================================
PROJECT_ROOT = os.environ.get("PROJECT_ROOT", "/path/to/GGN_progression_data")
PROJ = PROJECT_ROOT.rstrip("/") + "/"

# Input 1: Table S1 (three sheets, 999 nodules in total)
S1 = PROJ + "manuscript/Table/Table S1.xlsx"

# Input 2: bridge table (the source of the pipeline's 999 rows, row-aligned with kinetics.csv)
BRIDGE = (PROJ + "omics/supplementary_clinical_analysis/figure_assembly/"
                 "rerun_999/input/TableS1_curated_binary_nodule_type_999.xlsx")

# Input 3: raw pipeline kinetics (45 columns; supplies the 41 columns Table S1 lacks)
KIN = (PROJ + "omics/supplementary_clinical_analysis/figure_assembly/"
              "rerun_999/code/kinetics.csv")

OUT_DIR = os.path.dirname(os.path.abspath(__file__)) + "/"
OUT = OUT_DIR + "kinetics.csv"

# Table S1 column name -> kinetics column name
MAP = [
    ("3D VDT (days)",        "VDT"),
    ("3D MDT (days)",        "MDT3d_total"),
    ("3D MDT (solid, days)", "MDT3d_solid"),
    ("Annualized 3D CTR",    "ann_CTR3d"),
]
MC = [k for _, k in MAP]

LOG = []


def log(s=""):
    print(s)
    LOG.append(str(s))


def norm(x):
    """Excel cells arrive as str: '4.0' -> '4'; blank / NaN -> ''"""
    s = str(x).strip()
    if s in ("nan", "None", "NaT", "<NA>"):
        return ""
    return s[:-2] if s.endswith(".0") else s


def fnum(x):
    """To float; blank / non-numeric / NaN -> None"""
    try:
        v = float(str(x).strip())
        return None if v != v else v
    except (TypeError, ValueError):
        return None


def clean_group(s):
    """The plotting pipeline carries a historical typo, GK253GO-growth"""
    s = str(s).strip()
    return "GGO-growth" if s == "GK253GO-growth" else s


# ======================================================================
# 1. Read the three inputs
# ======================================================================
sheets = pd.read_excel(S1, sheet_name=None, dtype=str)
T = pd.concat([d.assign(__sheet=s) for s, d in sheets.items()], ignore_index=True)
log("Table S1 %s" % os.path.basename(S1))
for s, d in sheets.items():
    log("   sheet %-16s %4d rows" % (s, len(d)))
log("   merged %d rows" % len(T))

B = pd.read_excel(BRIDGE, dtype=str)
K = list(csv.DictReader(open(KIN, encoding="utf-8-sig")))
COLS = list(K[0].keys())
log("bridge table %d rows ; pipeline kinetics %d rows x %d columns" % (len(B), len(K), len(COLS)))
assert len(T) == len(B) == len(K) == 999, "all three sources must have 999 rows"


# ======================================================================
# 2. Two-hop matching: Table S1 --fingerprint--> bridge table --row index--> kinetics
#    Fingerprint = (follow-up days, baseline type, final type, pathology, group).
#    These five fields are encoded identically in both tables and are the only
#    reliable nodule-level key in this project.
# ======================================================================
def key_S1(t):
    r = T.loc[t]
    return (fnum(r["Follow-up duration (days)"]),
            norm(r["Baseline nodule type (0 = pGGO, 1 = PSN)"])[:1],
            norm(r["Final nodule type (0 = pGGO, 1 = PSN)"])[:1],
            norm(r["Pathological type (1 = AAH, 2 = AIS, 3 = MIA, 4 = IAC)"]).replace("/", ""),
            clean_group(r["Group"]))


def key_B(i):
    # The bridge table keeps its Chinese column headers, so these five fields are
    # read by their Chinese names; they hold the same values as key_S1 above.
    r = B.iloc[i]
    return (fnum(r["随访时长(天)"]),
            norm(r["基线结节性质"])[:1],
            norm(r["末次结节性质"])[:1],
            norm(r["病理类型"]).replace("/", ""),
            clean_group(r["Group"]))


pool = collections.defaultdict(list)
for i in range(len(B)):
    pool[key_B(i)].append(i)

pairs = {}
used = set()
for t in range(len(T)):
    free = [i for i in pool.get(key_S1(t), []) if i not in used]
    if free:
        pairs[t] = free[0]
        used.add(free[0])

log("")
log("Table S1 -> bridge table one-to-one pairing: %d / %d" % (len(pairs), len(T)))
if len(pairs) != len(T):
    miss = [t for t in range(len(T)) if t not in pairs]
    log("   [WARN] unmatched rows: %d : %s" % (len(miss), miss[:20]))
    raise SystemExit("incomplete pairing - stopping here rather than plotting a partial result")

# Spot check: the four kinetic columns should agree closely on both sides
same = tot = 0
for t, i in pairs.items():
    for sc, kc in MAP:
        a, b = fnum(T.loc[t, sc]), fnum(K[i].get(kc))
        if a is not None and b is not None:
            tot += 1
            same += abs(a - b) < 1e-6
log("   spot check: the four columns agree cell by cell %d / %d = %.2f%%" % (same, tot, 100.0 * same / tot))
log("   (the cells that disagree are exactly the nodules where the pipeline used a 2D back-filled value; they are overwritten by the Table S1 value below)")


# ======================================================================
# 3. Assemble: rows = Table S1 order; four kinetic columns + identity columns =
#    Table S1; the remaining 41 columns = pipeline
# ======================================================================
NEW = []
for t in range(len(T)):
    r = T.loc[t]
    src = K[pairs[t]]
    o = {c: src.get(c, "") for c in COLS}

    # -- identity / group / outcome: from Table S1
    o["mrn"] = norm(r["ID"])
    o["patient_id"] = norm(r["Patient ID"])
    o["group"] = clean_group(r["Group"])
    o["patho"] = norm(r["Pathological type (1 = AAH, 2 = AIS, 3 = MIA, 4 = IAC)"])
    o["outcome"] = norm(r["Clinical outcome (0 = surveillance, 1 = surgery)"])

    # -- follow-up: from Table S1
    fu = fnum(r["Follow-up duration (days)"])
    o["fu"] = "" if fu is None else norm(r["Follow-up duration (days)"])
    o["years"] = "" if not fu else repr(fu / 365.25)

    # -- the four kinetic columns: from Table S1. Blank stays blank; the pipeline
    #    value is **not** back-filled (back-filling would drag in the pipeline's 22
    #    stale 2D values and turn the pathology-subset VDT into 208, not 206)
    for sc, kc in MAP:
        v = fnum(r[sc])
        o[kc] = "" if v is None else repr(v)

    NEW.append(o)


# ======================================================================
# 4. Write out + reconciliation log
# ======================================================================
with open(OUT, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS)
    w.writeheader()
    w.writerows(NEW)

log("")
log("=" * 78)
log("write: %s  (%d rows x %d columns)" % (OUT, len(NEW), len(COLS)))

blank = lambda x: str(x).strip() in ("", "nan", "None", "/", "<NA>")
patho_set = [r for r in NEW if not blank(r.get("patho"))]
log("groups: " + "  ".join("%s=%d" % kv for kv in
                        sorted(collections.Counter(r["group"] for r in NEW).items())))
log("pathology non-missing = %d  (Table S1: 261)" % len(patho_set))
log("")
log("%-16s %-22s %-22s" % ("metric", "all: n/>0", "pathology: n/>0"))
for kc in MC:
    allv = [fnum(r[kc]) for r in NEW]
    allv = [x for x in allv if x is not None]
    pv = [fnum(r[kc]) for r in patho_set]
    pv = [x for x in pv if x is not None]
    log("%-16s %-22s %-22s" % (
        kc,
        "%4d / %4d" % (len(allv), sum(1 for x in allv if x > 0)),
        "%4d / %4d" % (len(pv), sum(1 for x in pv if x > 0))))

log("")
log("Reference: the pathology subset as computed from Table S1 itself (n=261) should be VDT 206 / MDT total 177 / MDT solid 112 / annualized CTR 134")
log("      whole table should be VDT 920/762, MDT total 920/583, MDT solid 346/255, annualized CTR 920/363")
log("")
log("Pipeline-only information retained:")
log("   three-state nodule type btype/ftype (ftype=2 means 'new solid component'; Table S1 has only two states)  ftype=2 total %d"
    % sum(1 for r in NEW if str(r.get("ftype")) == "2"))
log("   seven-state change pattern (including the hGGO intermediate state)  %s"
    % dict(collections.Counter(r["pattern"] for r in NEW)))

with open(OUT_DIR + "build_reconciliation_log.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(LOG) + "\n")
log("")
log("reconciliation log -> build_reconciliation_log.txt")
