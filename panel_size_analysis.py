import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt


def selectivity_entropy(profiles, baseline, epsilon=1e-10):
    shifted = np.maximum(profiles - baseline, 0)
    row_sums = np.where(shifted.sum(axis=1, keepdims=True) == 0, epsilon,
                        shifted.sum(axis=1, keepdims=True))
    p = shifted / row_sums
    return -(p * np.where(p > 0, np.log2(p + epsilon), 0)).sum(axis=1)


def gini_selectivity(profiles, baseline):
    shifted = np.maximum(profiles - baseline, 0)
    ginis = []
    for row in shifted:
        row_sorted = np.sort(row)
        n = len(row_sorted)
        total = row_sorted.sum()
        if total == 0:
            ginis.append(0.0); continue
        ginis.append((2 * np.sum((np.arange(1, n + 1) * row_sorted))) / (n * total) - (n + 1) / n)
    return np.array(ginis)


def s_score(profiles, threshold):
    return -(profiles > threshold).astype(float).mean(axis=1)


def ratio_selectivity(profiles, floor, top_n=1):
    ratios = []
    for row in profiles:
        s = np.sort(row)[::-1]
        ratios.append(s[0] - max(s[top_n] if len(s) > top_n else floor, floor))
    return np.array(ratios)


def to_ranks(scores):
    return len(scores) - scores.argsort().argsort()


def scores_at(M, baseline, threshold, floor):
    return {
        'entropy': to_ranks(-selectivity_entropy(M, baseline)),
        'gini':    to_ranks(gini_selectivity(M, baseline)),
        's_score': to_ranks(s_score(M, threshold)),
        'ratio':   to_ranks(ratio_selectivity(M, floor)),
    }


def panel_curves(M, baseline, threshold, floor, n_repeats=50, seed=42):
    """Mean/std Spearman r vs the full-panel ranking, over random subpanels."""
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


# Two sparse-coverage datasets (each compound profiled against a subset of the
# kinome), the appropriate substrate for a panel-size question. Davis and
# Anastassiadis are dense complete-matrix screens, so subsampling their kinases
# confounds panel size with apparent-active-count changes (see paper), and they
# are excluded here.
DATASETS = [
    dict(name='Klaeger', file='klaeger_matrix.csv', baseline=5.0, threshold=6.0, floor=5.0),
    dict(name='Metz',    file='metz_matrix.csv',    baseline=4.0, threshold=6.0, floor=4.0),
]

order = ['entropy', 'gini', 's_score', 'ratio']
colors = {'entropy': 'darkorange', 'gini': 'green', 's_score': 'steelblue', 'ratio': 'crimson'}

fig, axes = plt.subplots(1, len(DATASETS), figsize=(7 * len(DATASETS), 5), squeeze=False)
for ax, ds in zip(axes[0], DATASETS):
    M = pd.read_csv(ds['file'], index_col=0).values
    n_drugs, n_kinases = M.shape
    panel_sizes, results = panel_curves(M, ds['baseline'], ds['threshold'], ds['floor'])

    print(f"\n=== {ds['name']} ({n_drugs} drugs x {n_kinases} kinases) ===")
    print("Panel size stability (mean Spearman r vs full panel):")
    print(f"{'Panel size':>12s} {'entropy':>10s} {'gini':>10s} {'s_score':>10s} {'ratio':>10s}")
    for ps in panel_sizes:
        print(f"{ps:>12d}" + "".join(f"  {np.mean(results[d][ps]):+.3f}" for d in order))
    print("Minimum panel size for mean r > 0.90:")
    for d in order:
        print(f"  {d}: {pstar(panel_sizes, results, d)} kinases")

    for d in order:
        means = [np.mean(results[d][ps]) for ps in panel_sizes]
        stds = [np.std(results[d][ps]) for ps in panel_sizes]
        ax.plot(panel_sizes, means, color=colors[d], linewidth=2,
                label=f"{d} (p*={pstar(panel_sizes, results, d)})")
        ax.fill_between(panel_sizes, [m - s for m, s in zip(means, stds)],
                        [m + s for m, s in zip(means, stds)], alpha=0.15, color=colors[d])
    ax.axhline(0.90, color='black', linestyle='--', linewidth=1)
    ax.set_xlabel('Panel size (number of kinases)')
    ax.set_ylabel('Spearman r vs full panel ranking')
    ax.set_title(f"{ds['name']} ({n_drugs} cpd $\\times$ {n_kinases} kin)")
    ax.legend(fontsize=8)
    ax.set_ylim(-0.2, 1.01)

fig.suptitle('Rank stability as a function of kinase panel size (50 subsamples per size)',
             fontsize=13)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('panel_size_stability.png', dpi=150, bbox_inches='tight')
print("\nSaved panel_size_stability.png")
