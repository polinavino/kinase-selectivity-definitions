import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import pandas as pd
import numpy as np
from scipy.stats import spearmanr, kendalltau
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from itertools import product

# ── Load data ──────────────────────────────────────────────────────────────
affinity = pd.read_csv("davis_affinity.csv")
drugs    = pd.read_csv("davis_drugs.csv")
proteins = pd.read_csv("davis_proteins.csv")

matrix = affinity.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity')
matrix = matrix.fillna(5.0)
M = matrix.values          # shape (68, 433)
n_drugs, n_kinases = M.shape


# ── Selectivity definitions ─────────────────────────────────────────────────

def s_score(profiles, threshold, concentration=None):
    """
    Karaman S-score: fraction of kinases inhibited above a threshold.
    Lower score = more selective (hits fewer kinases).

    threshold: pKd value above which a kinase is considered 'hit'
    This is the original Karaman definition, parameterised by threshold.
    The threshold corresponds to a Kd cutoff — e.g. threshold=7 means Kd<100nM.
    """
    hits = (profiles > threshold).astype(float)
    return hits.mean(axis=1)  # fraction of kinases hit per drug


def selectivity_entropy(profiles, baseline=5.0, epsilon=1e-10):
    """
    Uitdehaag & Zaman selectivity entropy.
    Treats normalised activity as a probability distribution and computes
    Shannon entropy. Low entropy = activity concentrated on few targets = selective.

    baseline: the background pKd value (non-binder level) to subtract before normalising.
    This parameterisation matters: a higher baseline removes more 'noise' before
    computing the distribution, affecting which compounds look selective.
    """
    shifted = np.maximum(profiles - baseline, 0)
    row_sums = shifted.sum(axis=1, keepdims=True)
    # avoid division by zero for completely inactive compounds
    row_sums = np.where(row_sums == 0, epsilon, row_sums)
    p = shifted / row_sums
    # Shannon entropy: H = -sum(p * log(p)), with 0*log(0) = 0
    log_p = np.where(p > 0, np.log2(p + epsilon), 0)
    entropy = -(p * log_p).sum(axis=1)
    return entropy  # higher entropy = less selective


def gini_selectivity(profiles, baseline=5.0):
    """
    Gini coefficient of the activity distribution.
    Gini=1: all activity on one kinase (maximally selective).
    Gini=0: activity equally distributed (maximally promiscuous).

    baseline: subtract before computing Gini, same parameterisation issue as entropy.
    """
    shifted = np.maximum(profiles - baseline, 0)
    ginis = []
    for row in shifted:
        row_sorted = np.sort(row)
        n = len(row_sorted)
        cumsum = np.cumsum(row_sorted)
        total = cumsum[-1]
        if total == 0:
            ginis.append(0.0)
            continue
        gini = (2 * np.sum((np.arange(1, n+1) * row_sorted))) / (n * total) - (n+1)/n
        ginis.append(gini)
    return np.array(ginis)


def ratio_selectivity(profiles, top_n=1):
    """
    Ratio of primary target activity to nth-best off-target.
    Higher ratio = more selective.

    top_n: which off-target to use as denominator.
    top_n=1: ratio to best off-target (most stringent).
    top_n=2: ratio to second-best off-target, etc.
    This parameterisation matters enormously in practice.
    """
    ratios = []
    for row in profiles:
        sorted_desc = np.sort(row)[::-1]
        primary = sorted_desc[0]
        if len(sorted_desc) <= top_n:
            ratios.append(0.0)
            continue
        off_target = sorted_desc[top_n]
        if off_target <= 5.0:
            ratios.append(primary - 5.0)  # no meaningful off-target
        else:
            ratios.append(primary - off_target)
    return np.array(ratios)


# ── Sweep parameter spaces ──────────────────────────────────────────────────

# S-score thresholds (pKd): 5.5 to 8.0 in steps of 0.25
# covers range from weak binders (Kd~300nM) to strong binders (Kd~10nM)
s_thresholds = np.arange(5.5, 8.25, 0.25)

# Entropy baselines: 5.0 to 6.5
# how aggressively we filter out weak binders before computing distribution
ent_baselines = np.arange(5.0, 6.75, 0.25)

# Gini baselines: same range as entropy
gini_baselines = np.arange(5.0, 6.75, 0.25)

# Ratio top_n: 1 to 5
ratio_top_ns = list(range(1, 6))

print("Computing selectivity scores across parameter space...")

# For each drug, compute rank under every parameter combination
# Store as dict: {(definition, param): array of shape (n_drugs,)}
all_scores = {}

for t in s_thresholds:
    scores = s_score(M, threshold=t)
    # S-score: lower = more selective, so we negate for consistent ranking
    all_scores[('s_score', round(t,2))] = -scores

for b in ent_baselines:
    scores = selectivity_entropy(M, baseline=b)
    all_scores[('entropy', round(b,2))] = -scores  # negate: lower entropy = more selective

for b in gini_baselines:
    scores = gini_selectivity(M, baseline=b)
    all_scores[('gini', round(b,2))] = scores  # higher gini = more selective, no negation

for n in ratio_top_ns:
    scores = ratio_selectivity(M, top_n=n)
    all_scores[('ratio', n)] = scores  # higher ratio = more selective

# Convert to rankings (rank 1 = most selective)
all_ranks = {}
for key, scores in all_scores.items():
    # argsort twice gives rank
    ranks = len(scores) - scores.argsort().argsort()
    all_ranks[key] = ranks

print(f"Computed {len(all_scores)} score vectors")
print(f"Each covers {n_drugs} drugs")

# ── Rank stability analysis ─────────────────────────────────────────────────

keys = list(all_ranks.keys())
n_configs = len(keys)
rank_matrix = np.array([all_ranks[k] for k in keys])  # shape (n_configs, n_drugs)

# For each drug: std of its rank across all configurations
# High std = definition-sensitive (unstable)
# Low std = definition-robust (stable)
rank_std = rank_matrix.std(axis=0)
rank_mean = rank_matrix.mean(axis=0)

print("\nRank stability across all parameter configurations:")
print(f"Mean rank std per drug: {rank_std.mean():.2f}")
print(f"Min rank std (most stable): {rank_std.min():.2f} (Drug {rank_std.argmin()})")
print(f"Max rank std (most unstable): {rank_std.max():.2f} (Drug {rank_std.argmax()})")

# ── Within-definition stability ─────────────────────────────────────────────

print("\nWithin-definition rank stability (std of rank as parameter varies):")
for defn in ['s_score', 'entropy', 'gini', 'ratio']:
    defn_keys = [k for k in keys if k[0] == defn]
    defn_ranks = np.array([all_ranks[k] for k in defn_keys])
    drug_stds = defn_ranks.std(axis=0)
    print(f"  {defn:12s}: mean rank std = {drug_stds.mean():.2f}, "
          f"max = {drug_stds.max():.2f} (Drug {drug_stds.argmax()})")

# ── Cross-definition agreement ──────────────────────────────────────────────

print("\nSpearman correlations between definitions (using median rank per definition):")
median_ranks = {}
for defn in ['s_score', 'entropy', 'gini', 'ratio']:
    defn_keys = [k for k in keys if k[0] == defn]
    defn_ranks = np.array([all_ranks[k] for k in defn_keys])
    median_ranks[defn] = np.median(defn_ranks, axis=0)

defns = ['s_score', 'entropy', 'gini', 'ratio']
print(f"{'':12s}", end="")
for d in defns:
    print(f"{d:12s}", end="")
print()
for d1 in defns:
    print(f"{d1:12s}", end="")
    for d2 in defns:
        r, p = spearmanr(median_ranks[d1], median_ranks[d2])
        print(f"{r:12.3f}", end="")
    print()

# ── Save results ─────────────────────────────────────────────────────────────

results_df = pd.DataFrame({
    'drug_index': np.arange(n_drugs),
    'rank_mean': rank_mean,
    'rank_std': rank_std,
    'rank_cv': rank_std / (rank_mean + 1e-10),
})
for defn in ['s_score', 'entropy', 'gini', 'ratio']:
    results_df[f'median_rank_{defn}'] = median_ranks[defn]

results_df.to_csv("selectivity_results.csv", index=False)
print("\nSaved selectivity_results.csv")

# ── Plots ─────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: rank std per drug (overall instability)
ax = axes[0, 0]
sorted_idx = np.argsort(rank_std)
ax.bar(range(n_drugs), rank_std[sorted_idx], color='steelblue', alpha=0.7)
ax.set_xlabel('Drug (sorted by stability)')
ax.set_ylabel('Rank std across all configurations')
ax.set_title('Overall rank instability per drug\n(higher = more definition-sensitive)')

# Plot 2: rank trajectories for most/least stable drugs
ax = axes[0, 1]
most_stable = rank_std.argmin()
least_stable = rank_std.argmax()
mid_stable = np.argsort(rank_std)[n_drugs // 2]
for drug_idx, label, color in [
    (most_stable, f'Drug {most_stable} (most stable)', 'green'),
    (mid_stable, f'Drug {mid_stable} (median stability)', 'orange'),
    (least_stable, f'Drug {least_stable} (least stable)', 'red'),
]:
    ranks_for_drug = [all_ranks[k][drug_idx] for k in keys]
    ax.plot(ranks_for_drug, alpha=0.7, label=label, color=color)
ax.set_xlabel('Configuration index')
ax.set_ylabel('Rank (1 = most selective)')
ax.set_title('Rank trajectories across configurations')
ax.legend(fontsize=8)

# Plot 3: cross-definition correlation heatmap
ax = axes[1, 0]
corr_matrix = np.zeros((4, 4))
for i, d1 in enumerate(defns):
    for j, d2 in enumerate(defns):
        r, _ = spearmanr(median_ranks[d1], median_ranks[d2])
        corr_matrix[i, j] = r
im = ax.imshow(corr_matrix, vmin=-1, vmax=1, cmap='RdYlGn')
ax.set_xticks(range(4))
ax.set_yticks(range(4))
ax.set_xticklabels(defns, rotation=45)
ax.set_yticklabels(defns)
for i in range(4):
    for j in range(4):
        ax.text(j, i, f'{corr_matrix[i,j]:.2f}', ha='center', va='center', fontsize=10)
plt.colorbar(im, ax=ax)
ax.set_title('Spearman correlation between\nmedian rankings by definition')

# Plot 4: scatter of two most-disagreeing definitions
ax = axes[1, 1]
# find the pair with lowest correlation
min_corr = 1.0
min_pair = (defns[0], defns[1])
for i, d1 in enumerate(defns):
    for j, d2 in enumerate(defns):
        if i >= j:
            continue
        r, _ = spearmanr(median_ranks[d1], median_ranks[d2])
        if r < min_corr:
            min_corr = r
            min_pair = (d1, d2)

d1, d2 = min_pair
ax.scatter(median_ranks[d1], median_ranks[d2], alpha=0.7, edgecolors='k', linewidths=0.3)
for i in range(n_drugs):
    if abs(median_ranks[d1][i] - median_ranks[d2][i]) > 20:
        ax.annotate(f'Drug {i}', (median_ranks[d1][i], median_ranks[d2][i]),
                   fontsize=7, alpha=0.8)
ax.set_xlabel(f'Rank by {d1}')
ax.set_ylabel(f'Rank by {d2}')
ax.set_title(f'Most disagreeing definitions\n{d1} vs {d2} (r={min_corr:.2f})')

plt.tight_layout()
plt.savefig('selectivity_analysis.png', dpi=150, bbox_inches='tight')
print("Saved selectivity_analysis.png")
