import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

matrix = pd.read_csv("klaeger_matrix.csv", index_col=0)
M = matrix.values
n_drugs, n_kinases = M.shape

def selectivity_entropy(profiles, baseline=5.0, epsilon=1e-10):
    shifted = np.maximum(profiles - baseline, 0)
    row_sums = np.where(shifted.sum(axis=1, keepdims=True)==0, epsilon,
                        shifted.sum(axis=1, keepdims=True))
    p = shifted / row_sums
    return -(p * np.where(p>0, np.log2(p+epsilon), 0)).sum(axis=1)

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

def s_score(profiles, threshold=6.0):
    return -(profiles > threshold).astype(float).mean(axis=1)

def ratio_selectivity(profiles, top_n=1):
    ratios = []
    for row in profiles:
        s = np.sort(row)[::-1]
        ratios.append(s[0] - max(s[top_n] if len(s)>top_n else 5.0, 5.0))
    return np.array(ratios)

def to_ranks(scores):
    return len(scores) - scores.argsort().argsort()

ref_ranks = {
    'entropy': to_ranks(-selectivity_entropy(M)),
    'gini':    to_ranks(gini_selectivity(M)),
    's_score': to_ranks(s_score(M)),
    'ratio':   to_ranks(ratio_selectivity(M)),
}

panel_sizes = list(range(50, n_kinases+1, 30)) + [n_kinases]
n_repeats = 50
np.random.seed(42)
results = {d: {ps: [] for ps in panel_sizes} for d in ref_ranks}

for ps in panel_sizes:
    for _ in range(n_repeats):
        idx = np.random.choice(n_kinases, ps, replace=False)
        M_sub = M[:, idx]
        sub_ranks = {
            'entropy': to_ranks(-selectivity_entropy(M_sub)),
            'gini':    to_ranks(gini_selectivity(M_sub)),
            's_score': to_ranks(s_score(M_sub)),
            'ratio':   to_ranks(ratio_selectivity(M_sub)),
        }
        for d in results:
            r, _ = spearmanr(ref_ranks[d], sub_ranks[d])
            results[d][ps].append(r)

print("Panel size stability (mean Spearman r vs full panel):")
print(f"{'Panel size':>12s} {'entropy':>10s} {'gini':>10s} {'s_score':>10s} {'ratio':>10s}")
for ps in panel_sizes:
    row = f"{ps:>12d}"
    for d in ['entropy','gini','s_score','ratio']:
        row += f"  {np.mean(results[d][ps]):+.3f}"
    print(row)

print("\nMinimum panel size for mean r > 0.90:")
for d in ['entropy','gini','s_score','ratio']:
    for ps in panel_sizes:
        if np.mean(results[d][ps]) > 0.90:
            print(f"  {d}: {ps} kinases")
            break

fig, ax = plt.subplots(figsize=(8, 5))
colors = {'entropy':'darkorange','gini':'green','s_score':'steelblue','ratio':'crimson'}
for d in ['entropy','gini','s_score','ratio']:
    means = [np.mean(results[d][ps]) for ps in panel_sizes]
    stds  = [np.std(results[d][ps])  for ps in panel_sizes]
    ax.plot(panel_sizes, means, color=colors[d], label=d, linewidth=2)
    ax.fill_between(panel_sizes,
                    [m-s for m,s in zip(means,stds)],
                    [m+s for m,s in zip(means,stds)],
                    alpha=0.15, color=colors[d])
ax.axhline(0.90, color='black', linestyle='--', linewidth=1, label='r = 0.90')
ax.set_xlabel('Panel size (number of kinases)')
ax.set_ylabel('Spearman r vs full panel ranking')
ax.set_title('Rank stability as a function of kinase panel size\n(Klaeger dataset, 50 subsamples per panel size)')
ax.legend()
ax.set_ylim(-0.2, 1.01)
plt.tight_layout()
plt.savefig('panel_size_stability.png', dpi=150, bbox_inches='tight')
print("Saved panel_size_stability.png")
