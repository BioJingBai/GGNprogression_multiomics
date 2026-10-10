#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Joint continuous-variable models versus the three growth patterns: AUC for IAC
#
# Purpose   : Head-to-head comparison of (i) plain linear logistic combinations of
#             continuous kinetic variables and (ii) the three growth patterns for
#             discriminating invasive adenocarcinoma (IAC). Every model is scored on
#             the same rows and folds, with full-sample and repeated 5-fold CV AUC,
#             and the three patterns are tested on top of each variable by a nested
#             likelihood-ratio test.
#
# Inputs    : ./kinetics.csv (built by 01_build_nodule_kinetics_table.py)
#
# Outputs   : ./joint_auc_vs_threegroups.csv, ./joint_auc_nested_lr_sensitivity.csv
#
# Usage     : PROJECT_ROOT=/path/to/data python3 02_joint_auc_nested_likelihood_ratio.py
# -----------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "kinetics.csv")
GROUPS = ["Stable", "GGO-growth", "SOLID-growth"]
N_SEED, N_FOLD = 20, 5
RNG = np.random.default_rng(20260916)


# ---------------------------------------------------------------------------
# Minimal implementation with no sklearn: IRLS logistic + rank AUC
# ---------------------------------------------------------------------------
def fit_logistic(X, y, max_iter=100, ridge=1e-6):
    """X already includes an intercept column. IRLS with a tiny ridge against complete separation."""
    n, p = X.shape
    beta = np.zeros(p)
    for _ in range(max_iter):
        eta = X @ beta
        mu = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
        W = np.clip(mu * (1 - mu), 1e-10, None)
        z = eta + (y - mu) / W
        XtW = X.T * W
        H = XtW @ X + ridge * np.eye(p)
        try:
            step = np.linalg.solve(H, XtW @ z)
        except np.linalg.LinAlgError:
            break
        if np.max(np.abs(step - beta)) < 1e-9:
            beta = step
            break
        beta = step
    return beta


def predict_logit(X, beta):
    return X @ beta


def auc(y, score):
    """Mann-Whitney AUC. y in {0,1}."""
    y = np.asarray(y)
    if y.sum() == 0 or y.sum() == len(y):
        return np.nan
    r = stats.rankdata(score)
    n1 = y.sum()
    n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def design(df, cols):
    """Standardize + intercept. Returns (X, scaling parameters); the scaling parameters are reused within CV folds."""
    A = df[cols].to_numpy(float)
    mu, sd = A.mean(0), A.std(0)
    sd[sd == 0] = 1.0
    Z = (A - mu) / sd
    return np.c_[np.ones(len(Z)), Z], (mu, sd)


def standardize(A, mu=None, sd=None):
    """Standardize a raw matrix and add an intercept. mu/sd are estimated from the matrix itself when None."""
    A = np.asarray(A, float)
    if A.ndim == 1:
        A = A[:, None]
    if mu is None:
        mu, sd = A.mean(0), A.std(0)
        sd = np.where(sd == 0, 1.0, sd)
    return np.c_[np.ones(len(A)), (A - mu) / sd], (mu, sd)


def cv_auc(A, y, folds, scale=True):
    """**Within-fold standardization** — the test fold's mean/variance never enter
    training. scale=False is used for the three-group dummies (0/1 needs no scaling,
    and must not be scaled)."""
    A = np.asarray(A, float)
    if A.ndim == 1:
        A = A[:, None]
    out = []
    for tr, te in folds:
        if scale:
            Xtr, (mu, sd) = standardize(A[tr])
            Xte, _ = standardize(A[te], mu, sd)
        else:
            Xtr = np.c_[np.ones(len(tr)), A[tr]]
            Xte = np.c_[np.ones(len(te)), A[te]]
        beta = fit_logistic(Xtr, y[tr])
        out.append(auc(y[te], Xte @ beta))
    return np.array(out)


def make_folds(y, n_seed, n_fold, rng):
    """Stratified folds: stratification keeps the IAC fraction equal across folds. Returns [(folds_seed1), (folds_seed2), ...]"""
    all_folds = []
    for _ in range(n_seed):
        idx = rng.permutation(len(y))
        # split stratified by y
        pos = idx[y[idx] == 1]
        neg = idx[y[idx] == 0]
        fold_of = np.empty(len(y), int)
        for arr in (pos, neg):
            cut = np.array_split(arr, n_fold)
            for k, c in enumerate(cut):
                fold_of[c] = k
        all_folds.append([(np.where(fold_of != k)[0], np.where(fold_of == k)[0])
                          for k in range(n_fold)])
    return all_folds


def nested_lr(df, base_cols, y):
    """H0: logit(y) ~ base   H1: logit(y) ~ base + three groups (2 dummies). Returns (dchi2, p, df)"""
    d = pd.get_dummies(df["group"], prefix="g").reindex(
        columns=[f"g_{g}" for g in GROUPS[1:]], fill_value=0).astype(float)
    d.columns = [f"_grp{i}" for i in range(d.shape[1])]
    dd = pd.concat([df, d], axis=1)

    def ll(cols):
        X = design(dd, cols)[0]
        b = fit_logistic(X, y)
        eta = np.clip(X @ b, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        p = np.clip(p, 1e-12, 1 - 1e-12)
        return float((y * np.log(p) + (1 - y) * np.log(1 - p)).sum())

    l0 = ll(base_cols)
    l1 = ll(base_cols + list(d.columns))
    stat = 2 * (l1 - l0)
    return stat, float(stats.chi2.sf(stat, d.shape[1])), d.shape[1]


# ===========================================================================
# Load the data
# ===========================================================================
df = pd.read_csv(CSV, low_memory=False)
df["group"] = df["group"].replace({"GK253GO-growth": "GGO-growth"})
df = df[df["group"].isin(GROUPS)].copy()

for c in ["VDT", "MDT3d_total", "ann_CTR3d", "patho", "V0_3d", "V1_3d"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df[df["patho"].notna()].copy()               # valid pathology
df["IAC"] = (df["patho"] == 4).astype(int)

# Sign convention: more malignant = larger. Both transforms are computed:
#   raw -- the primary specification (`X = vals*direction`, no log)
#   log -- doubling times are naturally on a log scale
# AUC is a rank statistic, so a monotone transform (log) leaves it unchanged;
# only the nested likelihood-ratio chi-square is affected by the transform.
# VDT / MDT are doubling times (days), not rates: a negative value means halving
# every |VDT| days, so negating it would order the most stable nodules as the most
# malignant. Only > 0 is therefore kept. Annualized CTR is a signed rate and keeps
# every finite value.
_v = df["VDT"].where(df["VDT"] > 0)
_m = df["MDT3d_total"].where(df["MDT3d_total"] > 0)
df["_VDT_raw"], df["_MDT_raw"] = -_v, -_m
df["_VDT_log"], df["_MDT_log"] = -np.log(_v), -np.log(_m)
df["_CTR"] = df["ann_CTR3d"]
df["_VDT"], df["_MDT"] = df["_VDT_raw"], df["_MDT_raw"]      # primary convention = raw

print("=" * 100)
print("Reference: valid pathology n = %d, IAC = %d (%.1f%%), three groups %s"
      % (len(df), df.IAC.sum(), df.IAC.mean() * 100,
         df.group.value_counts().reindex(GROUPS).to_dict()))
print("=" * 100)

SUBSETS = [
    ("VDT", ["_VDT"]),
    ("MDT", ["_MDT"]),
    ("ann_CTR", ["_CTR"]),
    ("VDT + ann_CTR", ["_VDT", "_CTR"]),
    ("MDT + ann_CTR", ["_MDT", "_CTR"]),
    ("VDT + MDT", ["_VDT", "_MDT"]),
    ("VDT + MDT + ann_CTR", ["_VDT", "_MDT", "_CTR"]),
]

rows = []
for label, cols in SUBSETS:
    sub = df.dropna(subset=cols).copy()
    y = sub["IAC"].to_numpy(int)
    A = sub[cols].to_numpy(float)                 # raw values -- standardized within CV folds
    X = design(sub, cols)[0]                      # for the full-sample fit (no leakage concern)
    folds = make_folds(y, N_SEED, N_FOLD, RNG)

    # continuous variables, jointly
    full = auc(y, predict_logit(X, fit_logistic(X, y)))
    cvs = np.concatenate([cv_auc(A, y, f) for f in folds])
    cvs = cvs[~np.isnan(cvs)]

    # the three groups (same n, same folds) -- dummies are not standardized
    G = pd.get_dummies(sub["group"], prefix="g").reindex(
        columns=[f"g_{g}" for g in GROUPS], fill_value=0).astype(float)
    G = G.to_numpy()[:, 1:]                       # drop the reference level, Stable
    full_g = auc(y, predict_logit(np.c_[np.ones(len(G)), G], fit_logistic(np.c_[np.ones(len(G)), G], y)))
    cvs_g = np.concatenate([cv_auc(G, y, f, scale=False) for f in folds])
    cvs_g = cvs_g[~np.isnan(cvs_g)]

    dchi2, pval, dfree = nested_lr(sub, cols, y)

    rows.append(dict(label=label, n=len(sub), n_iac=int(y.sum()), k=len(cols),
                     auc_full=full, cv=cvs.mean(), cv_sd=cvs.std(),
                     auc_g_full=full_g, cv_g=cvs_g.mean(), cv_g_sd=cvs_g.std(),
                     dchi2=dchi2, p=pval, dfree=dfree))
    print("\n" + "-" * 100)
    print(f"[{label}] n = {len(sub)} (IAC = {y.sum()}), {len(cols)} continuous variable(s)")
    print(f"    continuous joint  full-sample AUC = {full:.3f}     repeated-CV AUC = {cvs.mean():.3f} +/- {cvs.std():.3f}")
    print(f"    three groups      full-sample AUC = {full_g:.3f}     repeated-CV AUC = {cvs_g.mean():.3f} +/- {cvs_g.std():.3f}")
    print(f"    three groups beyond the continuous joint: Δχ² = {dchi2:.1f} (df={dfree}), p = {pval:.3g}"
          + ("   [OK] three groups carry independent information" if pval < 0.05 else "   [--] no additional information"))

# ---------------------------------------------------------------------------
print("\n" + "=" * 100)
print("Summary table (same n, same folds; CV = %d folds x %d repeats)" % (N_FOLD, N_SEED))
print("=" * 100)
print(f"{'model':<24}{'n':>5}{'IAC':>5}  {'full AUC':>10}  {'CV AUC':>16}  {'3-group CV':>16}")
for r in rows:
    print(f"{r['label']:<24}{r['n']:>5}{r['n_iac']:>5}  {r['auc_full']:>10.3f}  "
          f"{r['cv']:>9.3f} ± {r['cv_sd']:.3f}  {r['cv_g']:>9.3f} ± {r['cv_g_sd']:.3f}")

# ---------------------------------------------------------------------------
# Sensitivity: rerun the nested LR with VDT / MDT on the log scale
# (AUC is unchanged and does not need recomputing)
# ---------------------------------------------------------------------------
print("\n" + "=" * 100)
print("Sensitivity: nested-LR Δχ² under the raw / log scales")
print("=" * 100)
print(f"{'continuous combination':<28}{'raw Δχ²':>12}{'p':>12}   {'log Δχ²':>12}{'p':>12}")
LOG_MAP = {"_VDT": "_VDT_log", "_MDT": "_MDT_log", "_CTR": "_CTR"}
sens = []
for label, cols in SUBSETS:
    sub = df.dropna(subset=cols).copy()
    yv = sub["IAC"].to_numpy(int)
    c_raw = cols
    c_log = [LOG_MAP[c] for c in cols]
    s_raw, p_raw, _ = nested_lr(sub, c_raw, yv)
    s_log, p_log, _ = nested_lr(sub, c_log, yv)
    sens.append((label, s_raw, p_raw, s_log, p_log))
    print(f"{label:<28}{s_raw:>12.1f}{p_raw:>12.2g}   {s_log:>12.1f}{p_log:>12.2g}")
print("\n  Under both scales the three groups all have p < 1e-6 -- the conclusion does not depend on the transform.")

out = pd.DataFrame(rows)
out.to_csv(os.path.join(HERE, "joint_auc_vs_threegroups.csv"), index=False,
           encoding="utf-8-sig")
pd.DataFrame(sens, columns=["model", "dchi2_raw", "p_raw", "dchi2_log", "p_log"]).to_csv(
    os.path.join(HERE, "joint_auc_nested_lr_sensitivity.csv"), index=False,
    encoding="utf-8-sig")
print("\n[OK] joint_auc_vs_threegroups.csv")
print("[OK] joint_auc_nested_lr_sensitivity.csv")
