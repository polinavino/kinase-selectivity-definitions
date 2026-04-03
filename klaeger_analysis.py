import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

matrix = pd.read_csv("klaeger_matrix.csv", index_col=0)
M = matrix.values
n_drugs, n_kinases = M.shape
drug_names = matrix.index.tolist()

def s_score(profiles, threshold):
    return -(profiles > threshold).astype(float).mean(axis=1)

def selectivity_entropy(profiles, baseline=5.0, epsilon=1e-10):
    shifted = np.maximum(profiles - baseline, 0)
    row_sums = np.where(shifted.sum(axis=1, keepdims=True)==0, epsilon, shifted.sum(axis=1, keepdims=True))
    p = shifted / row_sums
    return -(-(p * np.where(p>0, np.log2(p+epsilon), 0)).sum(axis=1))

def gini_selectivity(profiles, baseline=5.0):
    shifted = np.maximum(profiles - baseline, 0)
    ginis = []
    for row in shifted:
        row_sorted = np.sort(row)
        n = len(row_sorted)
        total = row_sorted.sum()
        if total == 0:
            ginis.append(0.0); continue
        ginis.append((2*np.sum((np.arange(1,n+1)*row_sorted)))/(n*total)-(n+1)/n)
    return np.array(ginis)

def ratio_selectivity(profiles, top_n=1):
    ratios = []
    for row in profiles:
        s = np.sort(row)[::-1]
        ratios.append(s[0] - max(s[top_n] if len(s)>top_n else 5.0, 5.0))
    return np.array(ratios)

def to_ranks(scores):
    return len(scores) - scores.argsort().argsort()

s_thresholds   = np.arange(5.5, 8.25, 0.25)
ent_baselines  = np.arange(5.0, 6.75, 0.25)
gini_baselines = np.arange(5.0, 6.75, 0.25)
ratio_top_ns   = list(range(1, 6))

s_ranks     = np.array([to_ranks(s_score(M, t))             for t in s_thresholds])
ent_ranks   = np.array([to_ranks(selectivity_entropy(M, b)) for b in ent_baselines])
gini_ranks  = np.array([to_ranks(gini_selectivity(M, b))    for b in gini_baselines])
ratio_ranks = np.array([to_ranks(ratio_selectivity(M, n))   for n in ratio_top_ns])

all_ranks  = np.vstack([s_ranks, ent_ranks, gini_ranks, ratio_ranks])
rank_std   = all_ranks.std(axis=0)
rank_mean  = all_ranks.mean(axis=0)

instability = {
    's_score': s_ranks.std(axis=0),
    'entropy': ent_ranks.std(axis=0),
    'gini':    gini_ranks.std(axis=0),
    'ratio':   ratio_ranks.std(axis=0),
}

median_ranks = {
    's_score': np.median(s_ranks, axis=0),
    'entropy': np.median(ent_ranks, axis=0),
    'gini':    np.median(gini_ranks, axis=0),
    'ratio':   np.median(ratio_ranks, axis=0),
}

print(f"Klaeger matrix: {n_drugs} drugs x {n_kinases} kinases")
print(f"\nRank stability:")
print(f"  Mean rank std: {rank_std.mean():.2f}")
print(f"  Most stable:   {drug_names[rank_std.argmin()]} ({rank_std.min():.2f})")
print(f"  Most unstable: {drug_names[rank_std.argmax()]} ({rank_std.max():.2f})")

print("\nWithin-definition stability:")
for name, ranks in [('s_score',s_ranks),('entropy',ent_ranks),('gini',gini_ranks),('ratio',ratio_ranks)]:
    stds = ranks.std(axis=0)
    print(f"  {name:12s}: mean={stds.mean():.2f}, max={stds.max():.2f} ({drug_names[stds.argmax()]})")

defns = ['s_score','entropy','gini','ratio']
print("\nSpearman correlations between definitions:")
print(f"{'':12s}", end="")
for d in defns: print(f"{d:12s}", end="")
print()
for d1 in defns:
    print(f"{d1:12s}", end="")
    for d2 in defns:
        r, _ = spearmanr(median_ranks[d1], median_ranks[d2])
        print(f"{r:12.3f}", end="")
    print()

n_active = (M > 6.0).sum(axis=1)
has_active = n_active > 0
print(f"\nZero-active drugs: {(~has_active).sum()}")
print(f"Rank std zero-active: {rank_std[~has_active].mean():.1f}")
print(f"Rank std active:      {rank_std[has_active].mean():.1f}")

features = {}
for i in range(n_drugs):
    profile = M[i]
    active = profile[profile > 6.0]
    s = np.sort(profile)[::-1]
    features[i] = {
        'n_active':      len(active),
        'top1_top2_gap': s[0]-s[1],
        'active_std':    active.std() if len(active)>1 else 0,
        'active_range':  (active.max()-active.min()) if len(active)>1 else 0,
        'top2_pkd':      s[1],
    }
feat_df = pd.DataFrame(features).T
mask = has_active

print("\nCorrelations with per-definition instability (active drugs only):")
print(f"{'Feature':20s} {'s_score':>10s} {'entropy':>10s} {'gini':>10s} {'ratio':>10s}")
for col in feat_df.columns:
    row = f"{col:20s}"
    for defn in defns:
        r, p = spearmanr(feat_df[col].values[mask], instability[defn][mask])
        sig = '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else ''
        row += f"  {r:+.3f}{sig:3s}"
    print(row)

disagreement = np.abs(median_ranks['ratio'] - median_ranks['entropy'])
print("\nTop 10 ratio-entropy disagreement drugs:")
print(f"{'Drug':22s} {'Disagree':>10s} {'n_active':>10s} {'top1_top2_gap':>15s}")
for idx in np.argsort(disagreement)[::-1][:10]:
    profile = M[idx]
    active = profile[profile>6.0]
    s = np.sort(profile)[::-1]
    print(f"{drug_names[idx]:22s} {disagreement[idx]:>10.1f} {len(active):>10d} {s[0]-s[1]:>15.3f}")

results_df = pd.DataFrame({
    'drug': drug_names,
    'rank_std': rank_std,
    'rank_mean': rank_mean,
    'median_rank_s_score': median_ranks['s_score'],
    'median_rank_entropy': median_ranks['entropy'],
    'median_rank_gini':    median_ranks['gini'],
    'median_rank_ratio':   median_ranks['ratio'],
})
results_df.to_csv("klaeger_selectivity_results.csv", index=False)
print("\nSaved klaeger_selectivity_results.csv")
