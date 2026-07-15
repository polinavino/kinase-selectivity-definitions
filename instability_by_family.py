"""
Figure 3: best predictor of within-family rank instability (Klaeger dataset).

For each definition family (S-score, entropy, Gini, ratio) we sweep the
definition's parameter, compute each active compound's rank standard deviation
across the sweep (its "instability"), correlate that instability against a set
of binding-profile features, and plot the single best-predicting feature per
family. Produces instability_by_family.png (Figure 3) and prints the full
feature x family correlation table (Table 2).
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

ACTIVE_THRESHOLD = 6.0   # pKd > 6 counts as an active kinase
BASELINE = 5.0           # assay detection floor

matrix = pd.read_csv("klaeger_matrix.csv", index_col=0)
M = matrix.values
n_drugs, n_kinases = M.shape
drug_names = matrix.index.tolist()


def s_score(profiles, threshold):
    return -(profiles > threshold).astype(float).mean(axis=1)


def selectivity_entropy(profiles, baseline=BASELINE, epsilon=1e-10):
    shifted = np.maximum(profiles - baseline, 0)
    row_sums = np.where(shifted.sum(axis=1, keepdims=True) == 0, epsilon,
                        shifted.sum(axis=1, keepdims=True))
    p = shifted / row_sums
    # Matches klaeger_analysis.py verbatim (sign convention affects tie-breaking
    # in the rank sweep, so reproduce it exactly to stay consistent with Table 2).
    return -(-(p * np.where(p > 0, np.log2(p + epsilon), 0)).sum(axis=1))


def gini_selectivity(profiles, baseline=BASELINE):
    shifted = np.maximum(profiles - baseline, 0)
    ginis = []
    for row in shifted:
        row_sorted = np.sort(row)
        n = len(row_sorted)
        total = row_sorted.sum()
        if total == 0:
            ginis.append(0.0)
            continue
        ginis.append((2 * np.sum(np.arange(1, n + 1) * row_sorted)) / (n * total) - (n + 1) / n)
    return np.array(ginis)


def ratio_selectivity(profiles, top_n=1):
    ratios = []
    for row in profiles:
        s = np.sort(row)[::-1]
        ratios.append(s[0] - max(s[top_n] if len(s) > top_n else BASELINE, BASELINE))
    return np.array(ratios)


def to_ranks(scores):
    return len(scores) - scores.argsort().argsort()


s_thresholds  = np.arange(5.5, 8.25, 0.25)
ent_baselines = np.arange(5.0, 6.75, 0.25)
gini_baselines = np.arange(5.0, 6.75, 0.25)
ratio_top_ns  = list(range(1, 6))

s_ranks     = np.array([to_ranks(s_score(M, t))             for t in s_thresholds])
ent_ranks   = np.array([to_ranks(selectivity_entropy(M, b)) for b in ent_baselines])
gini_ranks  = np.array([to_ranks(gini_selectivity(M, b))    for b in gini_baselines])
ratio_ranks = np.array([to_ranks(ratio_selectivity(M, n))   for n in ratio_top_ns])

instability = {
    's_score': s_ranks.std(axis=0),
    'entropy': ent_ranks.std(axis=0),
    'gini':    gini_ranks.std(axis=0),
    'ratio':   ratio_ranks.std(axis=0),
}

# Per-compound binding-profile features.
features = {}
for i in range(n_drugs):
    profile = M[i]
    active = profile[profile > ACTIVE_THRESHOLD]
    s = np.sort(profile)[::-1]
    features[i] = {
        'n_active':      len(active),
        'top1_top2_gap': s[0] - s[1],
        'active_std':    active.std() if len(active) > 1 else 0.0,
        'active_range':  (active.max() - active.min()) if len(active) > 1 else 0.0,
        'top2_pkd':      s[1],
    }
feat_df = pd.DataFrame(features).T
feature_cols = ['n_active', 'top1_top2_gap', 'active_std', 'active_range', 'top2_pkd']

# Restrict to active compounds (those with >=1 active kinase).
has_active = (M > ACTIVE_THRESHOLD).sum(axis=1) > 0
mask = has_active
n_active_drugs = int(mask.sum())

families = ['s_score', 'entropy', 'gini', 'ratio']
colors = {'s_score': 'tab:blue', 'entropy': 'tab:orange',
          'gini': 'tab:green', 'ratio': 'tab:red'}

# Full correlation table (Table 2) and best predictor per family.
print(f"Klaeger: {n_drugs} drugs x {n_kinases} kinases; active drugs n={n_active_drugs}\n")
print(f"{'Feature':16s} " + "".join(f"{f:>12s}" for f in families))
best = {}
corr = {f: {} for f in feature_cols}
for col in feature_cols:
    row = f"{col:16s} "
    for fam in families:
        r, p = spearmanr(feat_df[col].values[mask], instability[fam][mask])
        corr[col][fam] = (r, p)
        row += f"  {r:+.3f}{'*' if p < 0.05 else ' '} "
    print(row)
for fam in families:
    col = max(feature_cols, key=lambda c: abs(corr[c][fam][0]))
    r, p = corr[col][fam]
    best[fam] = (col, r, p)
    print(f"best predictor {fam:10s}: {col:14s} r={r:+.3f} p={p:.3f}")

# Plot: one panel per family, best-predicting feature vs. instability.
fig, axes = plt.subplots(1, 4, figsize=(20, 5))
for ax, fam in zip(axes, families):
    col, r, p = best[fam]
    x = feat_df[col].values[mask]
    y = instability[fam][mask]
    ax.scatter(x, y, color=colors[fam], edgecolors='white', linewidths=0.5, alpha=0.8, s=60)
    ax.set_title(f"{fam}\nbest feature: {col}\nr={r:.2f}, p={p:.3f}")
    ax.set_xlabel(col)
    ax.set_ylabel("Rank std")

fig.suptitle("Best predictor of instability per definition family", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("instability_by_family.png", dpi=150, bbox_inches='tight')
print("\nSaved instability_by_family.png")
