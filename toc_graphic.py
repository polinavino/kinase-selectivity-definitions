"""Graphical abstract / TOC graphic.

Rebuilds paper/figures/toc_graphic.png from current data. Left panel reuses
the Davis cross-definition correlation analysis (additional_datasets_analysis.py
load_davis/analyze) so it matches Table 3 exactly. Right panel repeats the gated
recovery analysis of Figure 3 (panel_recovery.py: same subpanel sizes, seed and
D1 gate) on Davis for three measures. Not one of the numbered manuscript
figures, so it is saved directly to paper/figures/ rather than via
script_logging's FIGURE_MAP.
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
    return -(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(axis=1)


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


# ------------------- right panel: gated recovery on Davis (Figure 6) -------------------
def neg_entropy(w):
    tot = w.sum(-1, keepdims=True)
    p = w / np.where(tot <= 0, 1.0, tot)
    return np.where(p > 0, p * np.log2(np.where(p > 0, p, 1.0)), 0.0).sum(-1)


X0, TAU = 5.0, 6.0
RIGHT = {  # label: (score, colour, marker), colours as in Figure 6
    "S-score, threshold at the floor": (lambda M: -(M > X0).mean(-1), "#6250d6", "X"),
    "entropy, published (1/Kd)": (lambda M: neg_entropy(np.where(M > X0, 10.0 ** (M - 9.0), 0.0)), "#e34948", "*"),
    "ratio R(1)": (lambda M: (lambda s: s[..., 0] - np.maximum(s[..., 1], X0))(-np.sort(-M, -1)), "#008300", "P"),
}
SIZES = [10, 20, 30, 40, 50, 65, 80, 100, 125, 150, 200, 250, 300]


def recovery_curves(M, n_repeats=50, seed=42):
    n, N = M.shape
    elig = (M > TAU).any(1)
    full = {k: f(M) for k, (f, _, _) in RIGHT.items()}
    sizes = [m for m in SIZES if m < N] + [N]
    rng = np.random.RandomState(seed)
    out = {k: [] for k in RIGHT}
    for m in sizes:
        acc = {k: [] for k in RIGHT}
        for _ in range(n_repeats if m < N else 1):
            Ms = M[:, rng.choice(N, m, replace=False)]
            both = (Ms > TAU).any(1) & elig
            for k, (f, _, _) in RIGHT.items():
                acc[k].append(spearmanr(full[k][both], f(Ms)[both])[0])
        for k in RIGHT:
            out[k].append(np.nanmean(acc[k]))
    return sizes, out


corr = davis_correlation_matrix()
aff = pd.read_csv(os.path.join(HERE, "davis_affinity.csv"))
Md = aff.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0).values
sizes, curves = recovery_curves(Md)
print(f"\nDavis gated recovery, smallest panel with mean Spearman > 0.90:")
for k in RIGHT:
    print(f"  {k}: {next((m for m, v in zip(sizes, curves[k]) if v > 0.90), None)} targets")

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
for k, (_, col, mk) in RIGHT.items():
    ax.plot(sizes, curves[k], color=col, marker=mk, ms=5, linewidth=2, label=k)
ax.axhline(0.90, color='#888888', linestyle='--', linewidth=1)
ax.set_xlabel('Panel size (number of targets)')
ax.set_ylabel('Spearman r with full-panel ranking')
ax.set_title('Recovery from random subpanels (Davis dataset)')
ax.legend(fontsize=8, loc='lower right', frameon=False)
ax.set_ylim(0, 1.01)

plt.tight_layout()
outpath = os.path.join(HERE, "paper", "figures", "toc_graphic.png")
plt.savefig(outpath, dpi=150, bbox_inches='tight')
print(f"\nSaved {outpath}")
