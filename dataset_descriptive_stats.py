"""
Dataset descriptive statistics cited in Methods but not previously saved to
any script output (flagged in the pre-resubmission audit): active fractions,
Klaeger zero-active-at-threshold and ratio-k-cap justification numbers,
Anastassiadis/Metz missing-entry percentages, and a robustness check for the
Davis stricter-cutoff claim (pKd >= 7 leaves cross-definition agreement
unchanged).
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))


def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')


def s_score(M, t):
    return -(M > t).astype(float).mean(axis=1)


def selectivity_entropy(M, b, eps=1e-10):
    sh = np.maximum(M - b, 0)
    rs = sh.sum(axis=1, keepdims=True)
    rs = np.where(rs == 0, eps, rs)
    p = sh / rs
    return -(p * np.where(p > 0, np.log2(p + eps), 0)).sum(axis=1)


def gini_selectivity(M, b):
    sh = np.maximum(M - b, 0)
    out = []
    for row in sh:
        rs = np.sort(row); n = len(rs); tot = rs.sum()
        out.append(0.0 if tot == 0 else (2 * np.sum(np.arange(1, n + 1) * rs)) / (n * tot) - (n + 1) / n)
    return np.array(out)


def ratio_selectivity(M, top_n, floor=5.0):
    out = []
    for row in M:
        s = np.sort(row)[::-1]
        off = s[top_n] if len(s) > top_n else floor
        out.append(s[0] - max(off, floor))
    return np.array(out)


DEFNS = ['s_score', 'entropy', 'gini', 'ratio']
BASELINES = np.arange(5.0, 6.75, .25)


def correlations(M, s_thr):
    ranks = {
        's_score': np.array([to_ranks(s_score(M, t)) for t in s_thr]),
        'entropy': np.array([to_ranks(-selectivity_entropy(M, b)) for b in BASELINES]),
        'gini':    np.array([to_ranks(gini_selectivity(M, b)) for b in BASELINES]),
        'ratio':   np.array([to_ranks(ratio_selectivity(M, k)) for k in range(1, 6)]),
    }
    med = {d: np.median(ranks[d], axis=0) for d in DEFNS}
    return np.array([[spearmanr(med[a], med[b])[0] for b in DEFNS] for a in DEFNS])


# --------------------------------- Davis -------------------------------------
aff = pd.read_csv(os.path.join(HERE, "davis_affinity.csv"))
M_davis = aff.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0).values
active_davis = (M_davis > 6.0).mean() * 100
print(f"Davis: active (pKd>6) fraction = {active_davis:.2f}% of {M_davis.size} entries")

corr_55_80 = correlations(M_davis, np.arange(5.5, 8.25, .25))
corr_65_90 = correlations(M_davis, np.arange(6.5, 9.25, .25))  # pKd>=7 stricter window, same width/step
print("\nDavis stricter-cutoff robustness check (S-score tau sweep shifted up by 1.0, "
      "from [5.5,8.0] to [6.5,9.0], i.e. a pKd>=7-centered window):")
print("Cross-definition correlation matrix, original tau sweep [5.5, 8.0]:")
for i, d in enumerate(DEFNS):
    print(f"  {d:8s}" + "".join(f"{corr_55_80[i, j]:8.3f}" for j in range(4)))
print("Cross-definition correlation matrix, stricter tau sweep [6.5, 9.0]:")
for i, d in enumerate(DEFNS):
    print(f"  {d:8s}" + "".join(f"{corr_65_90[i, j]:8.3f}" for j in range(4)))
print(f"Max absolute change in any off-diagonal cell: "
      f"{np.max(np.abs(corr_55_80 - corr_65_90)[~np.eye(4, dtype=bool)]):.3f}")

active7_davis = (M_davis > 7.0).mean() * 100
active6_davis = (M_davis > 6.0).mean() * 100
n_active_per_compound_6 = (M_davis > 6.0).sum(axis=1)
n_active_per_compound_7 = (M_davis > 7.0).sum(axis=1)
print(f"\nDavis active fraction at pKd>6: {active6_davis:.2f}%; at pKd>=7: {active7_davis:.2f}%")
print(f"Davis per-compound active-kinase count, mean at pKd>6: {n_active_per_compound_6.mean():.1f}; "
      f"at pKd>=7: {n_active_per_compound_7.mean():.1f} "
      f"(shrinkage factor {n_active_per_compound_6.mean() / n_active_per_compound_7.mean():.2f}x)")
print(f"Davis compounds pushed to zero active kinases by the stricter cutoff: "
      f"{(n_active_per_compound_7 == 0).sum()} of {M_davis.shape[0]} "
      f"(0 at pKd>6 cutoff: {(n_active_per_compound_6 == 0).sum()})")

# -------------------------------- Klaeger --------------------------------------
M_klaeger = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0).values
active_klaeger = (M_klaeger > 6.0).mean() * 100
print(f"\nKlaeger: active (pKd>6) fraction = {active_klaeger:.2f}% of {M_klaeger.size} entries")

n_active_6 = (M_klaeger > 6.0).sum(axis=1)
zero_active_6_pct = (n_active_6 == 0).mean() * 100
median_active = np.median(n_active_6)
pct_fewer_than_5 = (n_active_6 < 5).mean() * 100
print(f"Klaeger: median active kinases per compound (pKd>6) = {median_active:.1f}; "
      f"{pct_fewer_than_5:.1f}% of compounds have fewer than 5")

for tau in [6.0, 7.0, 8.0]:
    n_active = (M_klaeger > tau).sum(axis=1)
    pct_zero = (n_active == 0).mean() * 100
    print(f"Klaeger: {pct_zero:.1f}% of compounds have zero kinases above pKd = {tau:.1f}")

# --------------------------- Anastassiadis / Metz -----------------------------
raw = pd.read_excel(os.path.join(HERE, "anastassiadis.xls"), sheet_name='Sheet1', header=None)
pct_raw = raw.iloc[3:, 1:].apply(pd.to_numeric, errors='coerce').values
n_total = pct_raw.size
n_missing = np.isnan(pct_raw).sum()
print(f"\nAnastassiadis: {n_missing} of {n_total} entries missing in the raw matrix "
      f"({100 * n_missing / n_total:.2f}%)")

mz = pd.read_excel(os.path.join(HERE, "metz.xls"), sheet_name=0, header=0)
meta = ['Cmpd_ID', 'PUBCHEM_SID', 'Canonical_Smiles', 'External_Cmpd_ID', 'External_Source',
        'Cluster', 'ClusterSize', 'Cluster_MCSS', 'Molecular_Weight', 'ALogP',
        'Num_H_Acceptors', 'Num_H_Donors', 'tPSA', 'Promiscuity_1uM']
kin_cols = [c for c in mz.columns if c not in meta]
sub = mz[kin_cols].apply(pd.to_numeric, errors='coerce')
n_floor = (sub.values == 4.0).sum()
n_measured = (~np.isnan(sub.values)).sum()
print(f"Metz: {n_floor} of {n_measured} measured entries are exactly pKi = 4.0 "
      f"({100 * n_floor / n_measured:.3f}%)")
coverage_full = (~np.isnan(sub.values)).sum(axis=1)
print(f"Metz: median per-compound coverage across all {len(sub)} profiled compounds "
      f"(before the >=50-kinase restriction) = {np.median(coverage_full):.0f} of {len(kin_cols)} kinases")

print("\nDone.")
