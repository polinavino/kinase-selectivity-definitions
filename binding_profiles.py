"""
Figure 4: representative Davis binding profiles for the three instability types.

Selects three compounds directly from the rank-instability analysis of the Davis
dataset (same parameter sweep as selectivity_analysis.py):
  - most stable   = smallest overall rank std across all 30 configurations,
  - most unstable = largest overall rank std,
  - most ratio-sensitive = largest rank std within the ratio family,
and plots each compound's kinase affinities sorted in descending order.
Produces binding_profiles.png (Figure 4).
"""
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASELINE = 5.0
ACTIVE_THRESHOLD = 6.0

affinity = pd.read_csv("davis_affinity.csv")
matrix = affinity.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(BASELINE)
M = matrix.values
n_drugs, n_kinases = M.shape


# Definitions copied verbatim from selectivity_analysis.py so the instability
# ranking (and hence the three selected drugs) reproduces exactly.
def s_score(profiles, threshold):
    return (profiles > threshold).astype(float).mean(axis=1)


def selectivity_entropy(profiles, baseline=BASELINE, epsilon=1e-10):
    shifted = np.maximum(profiles - baseline, 0)
    row_sums = shifted.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, epsilon, row_sums)
    p = shifted / row_sums
    log_p = np.where(p > 0, np.log2(p + epsilon), 0)
    return -(p * log_p).sum(axis=1)


def gini_selectivity(profiles, baseline=BASELINE):
    shifted = np.maximum(profiles - baseline, 0)
    ginis = []
    for row in shifted:
        row_sorted = np.sort(row)
        n = len(row_sorted)
        total = np.cumsum(row_sorted)[-1]
        if total == 0:
            ginis.append(0.0)
            continue
        ginis.append((2 * np.sum(np.arange(1, n + 1) * row_sorted)) / (n * total) - (n + 1) / n)
    return np.array(ginis)


def ratio_selectivity(profiles, top_n=1):
    ratios = []
    for row in profiles:
        sorted_desc = np.sort(row)[::-1]
        primary = sorted_desc[0]
        if len(sorted_desc) <= top_n:
            ratios.append(0.0)
            continue
        off_target = sorted_desc[top_n]
        ratios.append(primary - 5.0 if off_target <= 5.0 else primary - off_target)
    return np.array(ratios)


def to_ranks(scores):
    return len(scores) - scores.argsort().argsort()


s_thresholds  = np.arange(5.5, 8.25, 0.25)
ent_baselines = np.arange(5.0, 6.75, 0.25)
gini_baselines = np.arange(5.0, 6.75, 0.25)
ratio_top_ns  = list(range(1, 6))

# Negate so rank 1 = most selective consistently across all four families,
# matching selectivity_analysis.py (otherwise the combined std is meaningless).
s_ranks     = np.array([to_ranks(-s_score(M, t))             for t in s_thresholds])
ent_ranks   = np.array([to_ranks(-selectivity_entropy(M, b)) for b in ent_baselines])
gini_ranks  = np.array([to_ranks(gini_selectivity(M, b))     for b in gini_baselines])
ratio_ranks = np.array([to_ranks(ratio_selectivity(M, n))    for n in ratio_top_ns])

all_ranks = np.vstack([s_ranks, ent_ranks, gini_ranks, ratio_ranks])
rank_std = all_ranks.std(axis=0)

most_stable    = int(rank_std.argmin())
most_unstable  = int(rank_std.argmax())
most_ratio_sen = int(ratio_ranks.std(axis=0).argmax())


def describe(idx):
    profile = M[idx]
    s = np.sort(profile)[::-1]
    n_active = int((profile > ACTIVE_THRESHOLD).sum())
    return profile, s, n_active


for label, idx in [("most stable", most_stable), ("most unstable", most_unstable),
                   ("most ratio-sensitive", most_ratio_sen)]:
    _, s, n_active = describe(idx)
    print(f"{label:22s}: Drug {idx}  sigma={rank_std[idx]:.2f}  "
          f"n_active={n_active}  top1-top2 gap={s[0]-s[1]:.3f}")

panels = [
    (most_stable,    "Most stable",          "clear winner", 'tab:green'),
    (most_unstable,  "Most unstable",        "flat profile", '#ef5350'),
    (most_ratio_sen, "Most ratio-sensitive", "tied top-2",   '#f0ad4e'),
]

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, (idx, title, note, color) in zip(axes, panels):
    profile, s, n_active = describe(idx)
    sorted_desc = np.sort(profile)[::-1]
    ax.bar(range(len(sorted_desc)), sorted_desc, color=color, width=1.0)
    ax.axhline(ACTIVE_THRESHOLD, ls='--', color='black', lw=1, label='pKd=6 threshold')
    ax.axhline(BASELINE, ls=':', color='gray', lw=1, label='baseline')
    ax.set_ylim(4.8, 11)
    ax.set_title(f"Drug {idx} — {title}\n({n_active} active kinases, {note})")
    ax.set_xlabel("Kinases (sorted by affinity)")
    ax.set_ylabel("pKd")
    ax.legend(loc='upper right')

plt.tight_layout()
plt.savefig("binding_profiles.png", dpi=150, bbox_inches='tight')
print("\nSaved binding_profiles.png")
