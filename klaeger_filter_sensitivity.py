"""
Klaeger high-confidence filter sensitivity (revision for Molecular Informatics,
Reviewer 1, major comment on filtered/imputed/undefined entries).

The manuscript matrix (klaeger_matrix.csv) keeps only the 5,123 drug-target
pairs annotated "High confidence" by Klaeger et al., out of 5,916 pairs with
any Apparent Kd measurement at all; the 793 "Low confidence" pairs are dropped
and, like every other untested pair, assigned the pKd=5.0 detection floor.
A true binding interaction removed by that filter is therefore folded into the
non-binding background rather than the active set.

This script reports (1) exactly how many entries the filter removes, and
(2) the sensitivity of the Table 1 Klaeger correlations and the panel-size p*
values to relaxing the filter (keeping High+Low confidence pairs instead of
High confidence only).
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


def ratio_selectivity(M, top_n, floor):
    out = []
    for row in M:
        s = np.sort(row)[::-1]
        off = s[top_n] if len(s) > top_n else floor
        out.append(s[0] - max(off, floor))
    return np.array(out)


DEFNS = ['s_score', 'entropy', 'gini', 'ratio']
S_THR = np.arange(5.5, 8.25, .25)
BASELINES = np.arange(5.0, 6.75, .25)


def correlations(M):
    ranks = {
        's_score': np.array([to_ranks(s_score(M, t)) for t in S_THR]),
        'entropy': np.array([to_ranks(-selectivity_entropy(M, b)) for b in BASELINES]),
        'gini':    np.array([to_ranks(gini_selectivity(M, b)) for b in BASELINES]),
        'ratio':   np.array([to_ranks(ratio_selectivity(M, k, 5.0)) for k in range(1, 6)]),
    }
    med = {d: np.median(ranks[d], axis=0) for d in DEFNS}
    return np.array([[spearmanr(med[a], med[b])[0] for b in DEFNS] for a in DEFNS])


def scores_at(M, baseline, threshold, floor):
    return {
        'entropy': to_ranks(-selectivity_entropy(M, baseline)),
        'gini':    to_ranks(gini_selectivity(M, baseline)),
        's_score': to_ranks(s_score(M, threshold)),
        'ratio':   to_ranks(ratio_selectivity(M, 1, floor)),
    }


def panel_curves(M, baseline, threshold, floor, n_repeats=50, seed=42):
    n_drugs, n_kinases = M.shape
    ref_ranks = scores_at(M, baseline, threshold, floor)
    panel_sizes = list(range(50, n_kinases, 30)) + [n_kinases]
    rng = np.random.RandomState(seed)
    results = {d: {ps: [] for ps in panel_sizes} for d in ref_ranks}
    for ps in panel_sizes:
        for _ in range(n_repeats):
            idx = rng.choice(n_kinases, ps, replace=False)
            sub = scores_at(M[:, idx], baseline, threshold, floor)
            for d in results:
                results[d][ps].append(spearmanr(ref_ranks[d], sub[d])[0])
    return panel_sizes, results


def pstar(panel_sizes, results, d):
    return next((ps for ps in panel_sizes if np.mean(results[d][ps]) > 0.90), None)


raw = pd.read_excel(os.path.join(HERE, "aan4368_Table_S2.xlsx"), sheet_name='Kinobeads')
n_total = len(raw)
n_high = (raw['Target Classification'] == 'High confidence').sum()
n_low = (raw['Target Classification'] == 'Low confidence').sum()
print(f"Raw Klaeger Kinobeads sheet: {n_total} drug-target pairs with an Apparent Kd")
print(f"  High confidence (kept in the manuscript matrix): {n_high}")
print(f"  Low confidence (dropped by the filter):          {n_low} "
      f"({100 * n_low / n_total:.1f}% of all measured pairs)")

raw['pKd'] = 9 - np.log10(raw['Apparent Kd'])
n_nan_kd = raw['pKd'].isna().sum()
if n_nan_kd:
    print(f"  Dropping {n_nan_kd} row(s) with a missing/non-numeric Apparent Kd "
          f"before aggregation (all {'High' if (raw.loc[raw['pKd'].isna(), 'Target Classification'] == 'High confidence').all() else 'Low'} confidence)")
    raw = raw.dropna(subset=['pKd'])

high_only = raw[raw['Target Classification'] == 'High confidence']
M_high = high_only.groupby(['Drug', 'Gene Name'])['pKd'].mean().unstack(fill_value=5.0)
M_relaxed = raw.groupby(['Drug', 'Gene Name'])['pKd'].mean().unstack(fill_value=5.0)
# Align to the same drug/kinase index (relaxed can only add, never remove, rows/cols)
M_relaxed = M_relaxed.reindex(index=M_high.index, columns=M_high.columns)
print(f"\nHigh-confidence-only matrix: {M_high.shape}")
print(f"Relaxed (High+Low confidence) matrix: {M_relaxed.shape}")
n_changed = int((np.abs(M_high.values - M_relaxed.values) > 1e-9).sum())
print(f"Cells changed by relaxing the filter: {n_changed} of {M_high.size} "
      f"({100 * n_changed / M_high.size:.2f}%)")

corr_high = correlations(M_high.values)
corr_relaxed = correlations(M_relaxed.values)
print("\nCross-definition Spearman correlation matrix (median ranks), High-confidence only:")
print(f"{'':8s}" + "".join(f"{d:8s}" for d in DEFNS))
for i, d in enumerate(DEFNS):
    print(f"{d:8s}" + "".join(f"{corr_high[i, j]:8.3f}" for j in range(4)))
print("\nSame, relaxed filter (High+Low confidence):")
print(f"{'':8s}" + "".join(f"{d:8s}" for d in DEFNS))
for i, d in enumerate(DEFNS):
    print(f"{d:8s}" + "".join(f"{corr_relaxed[i, j]:8.3f}" for j in range(4)))
print("\nMax absolute change in any off-diagonal cell: "
      f"{np.max(np.abs(corr_high - corr_relaxed)[~np.eye(4, dtype=bool)]):.3f}")

panel_sizes_h, results_h = panel_curves(M_high.values, baseline=5.0, threshold=6.0, floor=5.0)
panel_sizes_r, results_r = panel_curves(M_relaxed.values, baseline=5.0, threshold=6.0, floor=5.0)
print("\nPanel-size p* (smallest panel reaching mean rho > 0.90), High-confidence only vs. relaxed:")
for d in ['entropy', 'gini', 's_score', 'ratio']:
    print(f"  {d:8s}: {pstar(panel_sizes_h, results_h, d)} vs. {pstar(panel_sizes_r, results_r, d)}")
