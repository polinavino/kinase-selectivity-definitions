"""Graphical abstract / TOC graphic.

Rebuilds paper/figures/toc_graphic.png (5.0 cm wide, 600 dpi) from current data. Top panel reuses
the Davis cross-definition correlation analysis (additional_datasets_analysis.py
load_davis/analyze) so it matches Table 3 exactly. Bottom panel repeats the gated
recovery analysis of Figure 3 (panel_recovery.py: same subpanel sizes, seed and
D1 gate) on Davis for three measures: the S-score at the floor, the gated
entropy (T = tau* - x0) and the published entropy. Not one of the numbered manuscript
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


# ------------------- right panel: gated recovery on Davis (Figure 3) -------------------
def neg_entropy(w):
    tot = w.sum(-1, keepdims=True)
    p = w / np.where(tot <= 0, 1.0, tot)
    return np.where(p > 0, p * np.log2(np.where(p > 0, p, 1.0)), 0.0).sum(-1)


X0, TAU = 5.0, 6.0
# label: (score, colour, marker). The gated entropy and the published entropy keep
# their Figure 3 colours. Figure 3's violet for the S-score at the floor is too close
# to the gated entropy's blue at this size, so it takes dark green here.
RIGHT = {
    "S-score, threshold at the floor": (lambda M: -(M > X0).mean(-1), "#008300", "X"),
    "gated entropy": (lambda M: neg_entropy(np.where(M > X0, np.exp((M - X0) / (TAU - X0)), 0.0)), "#2a78d6", "o"),
    "entropy, published (1/Kd)": (lambda M: neg_entropy(np.where(M > X0, 10.0 ** (M - 9.0), 0.0)), "#e34948", "*"),
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

# Drawn at its print size: 5 cm wide (the journal's graphical-abstract width),
# so every font size below is the size on the page.
CM = 1 / 2.54
plt.rcParams.update({"font.size": 6, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "xtick.major.size": 2, "ytick.major.size": 2})
INK, MUTED = "#0b0b0b", "#52514e"
GRAYS = ["#ecebe8", "#c9c8c3", "#a3a29c", "#7d7c76", "#585752", "#3a3936"]
W, H = 5.0, 5.6   # cm; the width is fixed by the journal
fig = plt.figure(figsize=(W * CM, H * CM))
axh = fig.add_axes([0.25, 0.655, 0.72, 0.23])
axr = fig.add_axes([0.25, 0.14, 0.72, 0.32])

# Top: the six distinct pairs (lower triangle), log-affinity entropy and Gini.
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list("grays", GRAYS)
tri = np.where(np.tril(np.ones((4, 4), bool), -1), corr, np.nan)[1:, :3]
axh.imshow(tri, cmap=cmap, vmin=0.3, vmax=1.0, aspect="auto")
for i in range(3):
    for j in range(3):
        if not np.isnan(tri[i, j]):
            axh.text(j, i, f"{tri[i, j]:.2f}", ha="center", va="center", fontsize=6,
                     color="white" if tri[i, j] > 0.8 else INK)
axh.set_xticks(range(3)); axh.set_xticklabels(["S-score", "Entropy", "Gini"])
axh.set_yticks(range(3)); axh.set_yticklabels(["Entropy", "Gini", "Ratio"])
axh.tick_params(length=0)
for s in axh.spines.values():
    s.set_visible(False)
fig.text(0.5, 0.99, "Rank agreement between scores\n(log-affinity forms, Davis)", ha="center", va="top", fontsize=6.5)

# Bottom: recovery over the compounds scorable on both panels (D1 applied).
NAME = {"S-score, threshold at the floor": "S-score at floor",
        "gated entropy": "gated entropy",
        "entropy, published (1/Kd)": "published entropy"}
for k, (_, col, mk) in RIGHT.items():
    axr.plot(sizes, curves[k], color=col, marker=mk, ms=2.5, lw=1, mew=0, label=NAME[k])
axr.legend(loc="lower right", fontsize=5.5, frameon=False, handlelength=1.0,
           handletextpad=0.4, borderaxespad=0.0, borderpad=0.1, labelspacing=0.2)
axr.axhline(0.90, color="#888888", ls="--", lw=0.5)
axr.set_xlim(0, 370); axr.set_ylim(0, 1.05)
axr.set_xticks([0, 100, 200, 300]); axr.set_yticks([0, 0.5, 1])
axr.set_xlabel("Targets in random subpanel", labelpad=1)
axr.set_ylabel("Spearman with\nfull-panel ranking", labelpad=1)
for s in ("top", "right"):
    axr.spines[s].set_visible(False)
fig.text(0.5, 0.565, "Recovery from random subpanels\n(scorable compounds, Davis)", ha="center", va="top", fontsize=6.5)

outpath = os.path.join(HERE, "paper", "figures", "toc_graphic.png")
plt.savefig(outpath, dpi=600)
print(f"\nSaved {outpath}")
