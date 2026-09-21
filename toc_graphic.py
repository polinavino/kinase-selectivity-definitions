"""Graphical abstract / TOC graphic.

Rebuilds paper/figures/toc_graphic.png from current data. Left panel reuses
the Davis cross-definition correlation analysis (additional_datasets_analysis.py
load_davis/analyze) so it matches Table 1 exactly. Right panel reuses the
Klaeger panel-size stability analysis (panel_size_analysis.py) so it matches
Figure 5 exactly. Not one of the seven numbered manuscript figures, so it is
saved directly to paper/figures/ rather than via script_logging's FIGURE_MAP.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def s_score(M, threshold):
    return -(M > threshold).astype(float).mean(axis=1)


def selectivity_entropy(M, baseline, eps=1e-10):
    sh = np.maximum(M - baseline, 0)
    rs = sh.sum(axis=1, keepdims=True)
    rs = np.where(rs == 0, eps, rs)
    p = sh / rs
    return -(p * np.where(p > 0, np.log2(p + eps), 0)).sum(axis=1)


def gini_selectivity(M, baseline):
    sh = np.maximum(M - baseline, 0)
    out = []
    for row in sh:
        rs = np.sort(row); n = len(rs); tot = rs.sum()
        out.append(0.0 if tot == 0
                   else (2 * np.sum(np.arange(1, n + 1) * rs)) / (n * tot) - (n + 1) / n)
    return np.array(out)


def ratio_selectivity(M, top_n, floor):
    out = []
    for row in M:
        s = np.sort(row)[::-1]
        off = s[top_n] if len(s) > top_n else floor
        out.append(s[0] - max(off, floor))
    return np.array(out)


def to_ranks(scores):
    return rankdata(-np.asarray(scores), method='average')


DEFNS = ['s_score', 'entropy', 'gini', 'ratio']
LABELS = ['S', 'Ent', 'Gini', 'Ratio']


# --------------------------- left panel: Davis corr -------------------------
def davis_correlation_matrix():
    aff = pd.read_csv(os.path.join(HERE, "davis_affinity.csv"))
    M = aff.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0).values
    s_thr = np.arange(5.5, 8.25, .25)
    baselines = np.arange(5.0, 6.75, .25)
    ratio_floor = 5.0
    ranks = {
        's_score': np.array([to_ranks(s_score(M, t)) for t in s_thr]),
        'entropy': np.array([to_ranks(-selectivity_entropy(M, b)) for b in baselines]),
        'gini':    np.array([to_ranks(gini_selectivity(M, b)) for b in baselines]),
        'ratio':   np.array([to_ranks(ratio_selectivity(M, k, ratio_floor)) for k in range(1, 6)]),
    }
    med = {d: np.median(ranks[d], axis=0) for d in DEFNS}
    corr = np.array([[spearmanr(med[a], med[b])[0] for b in DEFNS] for a in DEFNS])
    print("Davis cross-definition Spearman correlations (median ranks):")
    print(f"{'':8s}" + "".join(f"{l:8s}" for l in LABELS))
    for i, d in enumerate(DEFNS):
        print(f"{LABELS[i]:8s}" + "".join(f"{corr[i, j]:8.3f}" for j in range(4)))
    return corr


# ------------------------ right panel: Klaeger stability --------------------
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


order = ['entropy', 'gini', 's_score', 'ratio']
colors = {'entropy': 'darkorange', 'gini': 'green', 's_score': 'steelblue', 'ratio': 'crimson'}

corr = davis_correlation_matrix()

km = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0)
M = km.values
n_drugs, n_kinases = M.shape
panel_sizes, results = panel_curves(M, baseline=5.0, threshold=6.0, floor=5.0)
print(f"\nKlaeger ({n_drugs} drugs x {n_kinases} kinases) minimum panel size for mean r > 0.90:")
for d in order:
    print(f"  {d}: {pstar(panel_sizes, results, d)} kinases")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
im = ax.imshow(corr, cmap='RdYlGn', vmin=0, vmax=1)
ax.set_xticks(range(4)); ax.set_xticklabels(LABELS)
ax.set_yticks(range(4)); ax.set_yticklabels(LABELS)
for i in range(4):
    for j in range(4):
        ax.text(j, i, f"{corr[i, j]:.2f}", ha='center', va='center', fontsize=10)
ax.set_title('Cross-definition agreement (Davis dataset)')
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

ax = axes[1]
for d in order:
    means = [np.mean(results[d][ps]) for ps in panel_sizes]
    ax.plot(panel_sizes, means, color=colors[d], linewidth=2, label=d)
ax.axhline(0.90, color='black', linestyle='--', linewidth=1)
ax.set_xlabel('Panel size (number of kinases)')
ax.set_ylabel('Spearman r vs full panel ranking')
ax.set_title('Rank stability vs panel size (Klaeger dataset)')
ax.legend(fontsize=8)
ax.set_ylim(-0.2, 1.01)

plt.tight_layout()
outpath = os.path.join(HERE, "paper", "figures", "toc_graphic.png")
plt.savefig(outpath, dpi=150, bbox_inches='tight')
print(f"\nSaved {outpath}")
