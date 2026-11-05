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
    return -(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(axis=1)


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


def candidate(M, floor, T=1.0, eps=1e-10):
    # gated smooth-hinge effective-target number (candidate_measure.py), selectivity orientation
    w = (M > floor) * (T * np.logaddexp(0.0, (M - floor) / T))
    rs = np.where(w.sum(1, keepdims=True) == 0, eps, w.sum(1, keepdims=True))
    p = w / rs
    return (p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(1)


def scores_at(M, baseline, threshold, floor):
    return {
        'entropy':   to_ranks(-selectivity_entropy(M, baseline)),
        'gini':      to_ranks(gini_selectivity(M, baseline)),
        's_score':   to_ranks(s_score(M, threshold)),
        'ratio':     to_ranks(ratio_selectivity(M, 1, floor)),
        'candidate': to_ranks(candidate(M, floor)),
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
M_full = raw.groupby(['Drug', 'Gene Name'])['pKd'].mean().unstack(fill_value=5.0)
# Analyses use the curated protein-kinase panel (curate_panels.py)
curated = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0).columns
print(f"Restricting every matrix to the {len(curated)} curated protein-kinase targets")
low_all = raw[raw['Target Classification'] != 'High confidence']
print(f"  low-confidence pairs on curated targets: {int(low_all['Gene Name'].isin(curated).sum())} of {len(low_all)}")
_lc = low_all[low_all['Gene Name'].isin(curated)]
print(f"  of these, fitted pKd above 5: {int((_lc['pKd'] > 5).sum())}; above 6: {int((_lc['pKd'] > 6).sum())}")
raw = raw[raw['Gene Name'].isin(curated)]
# measured values below the floor are set to the floor, as in curate_panels.py
M_high = M_high[curated].clip(lower=5.0)
M_full = M_full[[c for c in curated if c in M_full.columns]].clip(lower=5.0)
# Main variant: relax the filter on the same drug/kinase grid, so only the filter
# changes, not the panel. Low-confidence pairs for drugs or kinases with no
# high-confidence measurement fall outside this grid.
M_relaxed = M_full.reindex(index=M_high.index, columns=M_high.columns)
low = raw[raw['Target Classification'] != 'High confidence']
low_pairs = low.groupby(['Drug', 'Gene Name']).size().index
in_grid = [(d, g) for d, g in low_pairs if d in M_high.index and g in M_high.columns]
print(f"\nHigh-confidence-only matrix: {M_high.shape}")
print(f"Relaxed (High+Low confidence) matrix on the same grid: {M_relaxed.shape}")
print(f"Low-confidence pairs inside the high-confidence grid: {len(in_grid)} of {len(low_pairs)} "
      f"(the rest involve {len(set(M_full.index) - set(M_high.index))} drugs and "
      f"{len(set(M_full.columns) - set(M_high.columns))} kinases with no high-confidence measurement)")
print(f"Full relaxed matrix, all low-confidence pairs retained: {M_full.shape}")
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
act_h = (M_high.values > 6.0).any(1)
act_r = (M_relaxed.values > 6.0).any(1)
both_act = act_h & act_r
ch_a, cr_a = correlations(M_high.values[both_act]), correlations(M_relaxed.values[both_act])
print(f"\nCompounds with an active kinase under both filters ({both_act.sum()}), re-ranked among themselves: "
      f"max absolute change in any off-diagonal cell: "
      f"{np.max(np.abs(ch_a - cr_a)[~np.eye(4, dtype=bool)]):.3f}")
corr_full = correlations(M_full.values)
print("\nSame, all low-confidence pairs retained (adds the drugs with no high-confidence measurement):")
print(f"{'':8s}" + "".join(f"{d:8s}" for d in DEFNS))
for i, d in enumerate(DEFNS):
    print(f"{d:8s}" + "".join(f"{corr_full[i, j]:8.3f}" for j in range(4)))
print("Max absolute change vs high-confidence only: "
      f"{np.max(np.abs(corr_high - corr_full)[~np.eye(4, dtype=bool)]):.3f}")

panel_sizes_h, results_h = panel_curves(M_high.values, baseline=5.0, threshold=6.0, floor=5.0)
panel_sizes_r, results_r = panel_curves(M_relaxed.values, baseline=5.0, threshold=6.0, floor=5.0)
print("\nPanel-size p* without the D1 gate (SI Section S5 design), High-confidence only vs. relaxed:")
for d in ['entropy', 'gini', 's_score', 'ratio', 'candidate']:
    print(f"  {d:9s}: {pstar(panel_sizes_h, results_h, d)} vs. {pstar(panel_sizes_r, results_r, d)}")


# Gated recovery (main text, Figure 6; panel_recovery.py): same subpanel sizes, seed and D1
# gate, Spearman over compounds scored on the full panel and the subpanel.
def neg_entropy(w):
    tot = w.sum(-1, keepdims=True)
    p = w / np.where(tot <= 0, 1.0, tot)
    return np.where(p > 0, p * np.log2(np.where(p > 0, p, 1.0)), 0.0).sum(-1)


GATED = {
    "S-score tau*": lambda M: -(M > 6.0).mean(-1),
    "S-score x0":   lambda M: -(M > 5.0).mean(-1),
    "entropy x0":   lambda M: neg_entropy(np.maximum(M - 5.0, 0.0)),
    "entropy tau*": lambda M: neg_entropy(np.maximum(M - 6.0, 0.0)),
    "Gini":         lambda M: gini_selectivity(M, 5.0),
    "ratio":        lambda M: (lambda s: s[..., 0] - np.maximum(s[..., 1], 5.0))(-np.sort(-M, -1)),
    "candidate":    lambda M: neg_entropy(np.where(M > 5.0, np.logaddexp(0.0, M - 5.0), 0.0)),
    "entropy Ka":   lambda M: neg_entropy(np.where(M > 5.0, 10.0 ** (M - 9.0), 0.0)),
    "entropy Ka, floor kept": lambda M: neg_entropy(10.0 ** (M - 9.0)),
    "partition Ka": lambda M: (lambda K: K.max(-1) / np.where(K.sum(-1) > 0, K.sum(-1), 1.0))(
        np.where(M > 5.0, 10.0 ** (M - 9.0), 0.0)),
}


def gated_pstar(M, sizes=(10, 20, 30, 40, 50, 65, 80, 100, 125, 150, 200), R=50, seed=42):
    n, N = M.shape
    elig = (M > 6.0).any(1)
    full = {k: f(M) for k, f in GATED.items()}
    rng = np.random.RandomState(seed)
    ps, cov50 = {k: None for k in GATED}, None
    for m in list(sizes) + [N]:
        acc = {k: [] for k in GATED}
        cov = []
        for _ in range(R if m < N else 1):
            Ms = M[:, rng.choice(N, m, replace=False)]
            both = (Ms > 6.0).any(1) & elig
            cov.append(both.sum() / elig.sum())
            for k, f in GATED.items():
                acc[k].append(spearmanr(full[k][both], f(Ms)[both])[0])
        if m == 50:
            cov50 = np.mean(cov)
        for k in GATED:
            if ps[k] is None and np.nanmean(acc[k]) > 0.90:
                ps[k] = m
    return ps, cov50


ph, ch = gated_pstar(M_high.values)
pr, cr = gated_pstar(M_relaxed.values)
print("\nGated recovery p* (main text, Figure 3), High-confidence only vs. relaxed:")
for k in GATED:
    print(f"  {k:22s}: {ph[k]} vs. {pr[k]}")
print(f"  coverage at 50 targets: {100 * ch:.1f}% vs. {100 * cr:.1f}%")
